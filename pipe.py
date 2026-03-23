#!/usr/bin/env python3
"""
在已有 system.gro（DPPC/DOPC/CHOL/DPPE/W/离子）基础上：
1) 在 DPPE/DSPE 的 NH3 上沿「NH3→球心」的外向法线生长 PEG；
2) 用 MDAnalysis 移除全部水与溶剂离子，并统一用 MDAnalysis 读写 gro，保证格式正确。
配置见 config.yaml。
"""
from __future__ import annotations

import math
import shutil
import subprocess
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None

try:
    import MDAnalysis as mda
    import numpy as np
    _has_mda = True
except ImportError:
    mda = None
    np = None
    _has_mda = False

VERSION3_DIR = Path(__file__).resolve().parent
VERSION2_DIR = VERSION3_DIR  # 扁平化后 ndx.py 与 pipe.py 同目录


def _resolve_path(p: str, base: Path) -> Path:
    if not p or not str(p).strip():
        return base
    path = Path(str(p).strip().strip('"').strip("'"))
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def _norm3(x, y, z):
    n = math.sqrt(x * x + y * y + z * z)
    if n < 1e-10:
        return 1.0, 0.0, 0.0
    return x / n, y / n, z / n


def _parse_gro(path: Path):
    """解析 GRO 文件，返回 title, residues, box。residues 为列表，每项为 (resnr, resname, atoms)，atoms 为 [(atname, x, y, z), ...]。坐标 nm。"""
    with open(path, "r", encoding="utf-8") as f:
        lines = [l for l in f]
    if len(lines) < 3:
        raise ValueError(f"GRO 行数过少: {path}")
    title = lines[0].rstrip()
    try:
        n_atoms = int(lines[1].strip())
    except ValueError:
        raise ValueError(f"GRO 第二行应为原子数: {path}")
    # 原子行：%5d%-5s%5s%5d%8.3f%8.3f%8.3f
    residues = []
    current = None
    for i in range(2, 2 + n_atoms):
        if i >= len(lines):
            break
        line = lines[i]
        if len(line) < 44:
            continue
        resnr = int(line[0:5])
        resname = line[5:10].strip()
        atname = line[10:15].strip()
        atid = int(line[15:20])
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        if current is None or (current[0], current[1]) != (resnr, resname):
            current = (resnr, resname, [])
            residues.append(current)
        current[2].append((atname, x, y, z))
    # box
    box = (40.0, 40.0, 40.0)
    if len(lines) > 2 + n_atoms:
        parts = lines[2 + n_atoms].split()
        if len(parts) >= 3:
            box = (float(parts[0]), float(parts[1]), float(parts[2]))
    return title, residues, box, n_atoms


def _sphere_center_from_config(config: dict, box) -> tuple[float, float, float]:
    center = config.get("sphere_center")
    if center is None or (isinstance(center, str) and center.strip().lower() == "auto"):
        return box[0] / 2.0, box[1] / 2.0, box[2] / 2.0
    if isinstance(center, (list, tuple)) and len(center) >= 3:
        return float(center[0]), float(center[1]), float(center[2])
    return box[0] / 2.0, box[1] / 2.0, box[2] / 2.0


def _which_dppe_to_peg(config: dict, dppe_resids: list[int]) -> set[int]:
    to_peg = config.get("dppe_to_peg", config.get("peg_resids", "all"))
    if to_peg == "all" or (isinstance(to_peg, str) and to_peg.strip().lower() == "all"):
        return set(dppe_resids)
    if isinstance(to_peg, (list, tuple)):
        return set(int(r) for r in to_peg)
    return set()


def _grow_peg_build_new_atoms(u, config, head_resname, convert_set, peg_resname, peg_length, bond, cx, cy, cz):
    """从 MDAnalysis Universe 构建 PEG 生长后的 (resid, resname, atname, x_nm, y_nm, z_nm) 列表。"""
    new_atoms = []
    for res in u.residues:
        resnr = res.resid
        resname = (res.resname or "").strip()
        resname_upper = resname.upper()
        if resname_upper == head_resname and resnr in convert_set:
            nh3_atom = None
            for at in res.atoms:
                if (at.name or "").strip().upper() == "NH3":
                    nh3_atom = at
                    break
            if nh3_atom is None:
                for at in res.atoms:
                    pos = at.position / 10.0  # Angstrom -> nm
                    new_atoms.append((resnr, peg_resname, at.name, pos[0], pos[1], pos[2]))
                continue
            nh3_pos = nh3_atom.position / 10.0  # nm
            nx, ny, nz = nh3_pos[0] - cx, nh3_pos[1] - cy, nh3_pos[2] - cz
            ux, uy, uz = _norm3(nx, ny, nz)
            for at in res.atoms:
                pos = at.position / 10.0
                new_atoms.append((resnr, peg_resname, at.name, pos[0], pos[1], pos[2]))
            for k in range(1, peg_length + 1):
                d = bond * k
                new_atoms.append((
                    resnr, peg_resname, "EC",
                    nh3_pos[0] + d * ux, nh3_pos[1] + d * uy, nh3_pos[2] + d * uz,
                ))
        else:
            for at in res.atoms:
                pos = at.position / 10.0
                new_atoms.append((resnr, resname, at.name, pos[0], pos[1], pos[2]))
    return new_atoms


def _write_gro_manual(new_atoms: list, box_nm: tuple, path: Path, title: str = "LNB with PEG") -> None:
    """按 GRO 格式将 (resid, resname, atname, x, y, z) 列表写入文件（用于后续由 MDAnalysis 再次写出）。
    GRO 规定 resnr 与 atom id 仅 5 位（1～99999），resnr 按残基顺序输出为 1,2,3,… 避免列错位。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    total = len(new_atoms)
    lines = [title + "\n", f"{total:5d}\n"]
    resnr_prev = None
    resnr_out = 0
    for i, (resnr, rname, aname, x, y, z) in enumerate(new_atoms, start=1):
        atid_5 = (i - 1) % 99999 + 1
        if resnr != resnr_prev:
            resnr_prev = resnr
            resnr_out = (resnr_out % 99999) + 1
        rname_5 = ((rname or "").strip()[:5]).ljust(5)
        aname_5 = ((aname or "").strip()[:5]).ljust(5)
        lines.append(f"{resnr_out:5d}{rname_5}{aname_5}{atid_5:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
    lines.append(f"{box_nm[0]:10.5f}{box_nm[1]:10.5f}{box_nm[2]:10.5f}\n")
    with open(path, "w", encoding="utf-8") as f:
        f.writelines(lines)


def _write_gro_with_mda(new_atoms: list, box_nm: tuple, output_gro: Path, title: str = "LNB with PEG") -> None:
    """用 MDAnalysis 写出 GRO：先按格式写临时文件，再由 MDAnalysis 读取并写入目标，保证格式正确。"""
    import tempfile
    with tempfile.NamedTemporaryFile(mode="w", suffix=".gro", delete=False) as f:
        tmp_path = Path(f.name)
    try:
        _write_gro_manual(new_atoms, box_nm, tmp_path, title=title)
        u = mda.Universe(str(tmp_path))
        output_gro.parent.mkdir(parents=True, exist_ok=True)
        u.atoms.write(str(output_gro))
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def grow_peg_on_dppe(
    input_gro: Path,
    output_gro: Path,
    config: dict,
) -> int:
    """
    在 DPPE/DSPE 的 NH3 上沿球心→NH3 的外向法线生长 PEG，写出带 PEG 的 gro。
    使用 MDAnalysis 读写 gro，保证格式正确。返回新总原子数。
    """
    if not _has_mda:
        # 无 MDAnalysis 时回退到手工解析与写入
        title, residues, box, n_atoms = _parse_gro(input_gro)
        cx, cy, cz = _sphere_center_from_config(config, box)
        peg_length = int(config.get("peg_length", 15))
        peg_resname = (config.get("peg_resname") or "15PEG").strip()[:5]
        bond = float(config.get("peg_bond_length_nm", 0.36))
        head_resname = (config.get("head_resname") or "DPPE").strip().upper()
        dppe_resids = [r[0] for r in residues if (r[1] or "").strip().upper() == head_resname]
        convert_set = _which_dppe_to_peg(config, dppe_resids)
        if not convert_set:
            import shutil
            shutil.copy2(input_gro, output_gro)
            return n_atoms
        new_atoms = []
        for resnr, resname, atoms in residues:
            resname_stripped = (resname or "").strip().upper()
            if resname_stripped == head_resname and resnr in convert_set:
                nh3_pos = None
                for atname, x, y, z in atoms:
                    if (atname or "").strip().upper() == "NH3":
                        nh3_pos = (x, y, z)
                        break
                if nh3_pos is None:
                    for a in atoms:
                        new_atoms.append((resnr, peg_resname, a[0], a[1], a[2], a[3]))
                    continue
                nx, ny, nz = nh3_pos[0] - cx, nh3_pos[1] - cy, nh3_pos[2] - cz
                ux, uy, uz = _norm3(nx, ny, nz)
                for atname, x, y, z in atoms:
                    new_atoms.append((resnr, peg_resname, atname, x, y, z))
                for k in range(1, peg_length + 1):
                    d = bond * k
                    new_atoms.append((resnr, peg_resname, "EC",
                        nh3_pos[0] + d * ux, nh3_pos[1] + d * uy, nh3_pos[2] + d * uz))
            else:
                for atname, x, y, z in atoms:
                    new_atoms.append((resnr, resname, atname, x, y, z))
        total = len(new_atoms)
        out_lines = [title + "\n", f"{total:5d}\n"]
        for i, (resnr, rname, aname, x, y, z) in enumerate(new_atoms, start=1):
            atid_5 = (i - 1) % 99999 + 1
            out_lines.append(f"{resnr:5d}{rname:5s}{aname:5s}{atid_5:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
        out_lines.append(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n")
        output_gro.parent.mkdir(parents=True, exist_ok=True)
        with open(output_gro, "w", encoding="utf-8") as f:
            f.writelines(out_lines)
        print(f"   已在 {len(convert_set)} 个 {head_resname} 上生长 PEG（无 MDAnalysis，手工写入），总原子数 {total}")
        return total

    u = mda.Universe(str(input_gro))
    dims = u.dimensions
    if dims is not None and len(dims) >= 3:
        box_nm = (dims[0] / 10.0, dims[1] / 10.0, dims[2] / 10.0)
    else:
        box_nm = (40.0, 40.0, 40.0)
    cx, cy, cz = _sphere_center_from_config(config, box_nm)
    peg_length = int(config.get("peg_length", 15))
    peg_resname = (config.get("peg_resname") or "15PEG").strip()[:5]
    bond = float(config.get("peg_bond_length_nm", 0.36))
    head_resname = (config.get("head_resname") or "DPPE").strip().upper()
    dppe_resids = [res.resid for res in u.residues if (res.resname or "").strip().upper() == head_resname]
    convert_set = _which_dppe_to_peg(config, dppe_resids)
    if not convert_set:
        print(f"  未配置要转为 PEG 的 {head_resname}，直接复制 gro。")
        import shutil
        shutil.copy2(input_gro, output_gro)
        return len(u.atoms)

    new_atoms = _grow_peg_build_new_atoms(u, config, head_resname, convert_set, peg_resname, peg_length, bond, cx, cy, cz)
    total = len(new_atoms)
    _write_gro_with_mda(new_atoms, box_nm, output_gro, title="LNB with PEG")
    print(f"   已在 {len(convert_set)} 个 {head_resname} 上生长 PEG，写出 {output_gro}，总原子数 {total}")
    return total


def remove_all_solvent(
    gro_path: Path,
    config: dict,
    out_gro: Path,
    out_top: Path | None,
) -> None:
    """移除全部水和溶剂离子（W/NA/CL），仅保留脂质、PEG、气体等，用于真空模拟。用 MDAnalysis 写出 gro。"""
    if not _has_mda:
        print("   跳过移除溶剂：未安装 MDAnalysis。pip install MDAnalysis")
        import shutil
        if out_gro != gro_path:
            shutil.copy2(gro_path, out_gro)
        return

    peg_resname = (config.get("peg_resname") or "15PEG").strip()[:5]
    solvent_set = {"W", "NA", "CL"}

    u = mda.Universe(str(gro_path))
    non_solvent = [r for r in u.residues if (r.resname or "").strip().upper() not in solvent_set]
    n_removed = len(u.residues) - len(non_solvent)
    print(f"   已移除全部水与溶剂离子，共 {n_removed} 个溶剂残基")

    # 按 [ molecules ] 顺序排列残基，与 top 一致
    order_first = ("DPPC", "DOPC", "DOPS", "CHOL", "DPPE", "DSPE", peg_resname, "O2", "W", "NA", "CL")
    by_name = {}
    for r in non_solvent:
        name = (r.resname or "").strip()
        if name not in by_name:
            by_name[name] = []
        by_name[name].append(r)
    rest_names = sorted(k for k in by_name if k not in order_first)
    ordered_residues = []
    for name in order_first + tuple(rest_names):
        ordered_residues.extend(by_name.get(name, []))

    # 用 MDAnalysis Merge 按顺序合并残基的 AtomGroup，resid 重排为 1,2,3,…（GRO 仅支持 5 位）
    merged = mda.Merge(*[r.atoms for r in ordered_residues])
    merged.dimensions = u.dimensions
    resids_new = np.array([(i - 1) % 99999 + 1 for i in range(1, len(merged.residues) + 1)], dtype=np.int32)
    merged.residues.resids = resids_new
    out_gro.parent.mkdir(parents=True, exist_ok=True)
    merged.atoms.write(str(out_gro))

    if out_top and out_top.suffix == ".top":
        _write_simple_top(out_top, ordered_residues, peg_resname)


def copy_peg_itp(config: dict, out_dir: Path) -> None:
    """若配置了 dpeg_itp，复制到 out_dir/itps 并改为 moleculetype 与 peg_resname 一致（如 45PEG）。
    参考 version2/pipe.py 的 copy_15peg_itp。"""
    src = config.get("dpeg_itp")
    if not src or not str(src).strip():
        return
    src_path = _resolve_path(src.strip(), VERSION3_DIR)
    if not src_path.is_file():
        print(f"   未找到 PEG itp 源文件: {src_path}，跳过复制")
        return
    peg_resname = (config.get("peg_resname") or "15PEG").strip()[:5]
    itps_dst = out_dir / "itps"
    itps_dst.mkdir(parents=True, exist_ok=True)
    dst_path = itps_dst / f"{peg_resname}.itp"
    with open(src_path, "r", encoding="utf-8") as f:
        text = f.read()
    # 将 [ moleculetype ] 下一行的分子名改为 peg_resname（如 DPEG/15PEG -> 45PEG）
    lines = text.split("\n")
    for i, line in enumerate(lines):
        if "[ moleculetype ]" in line:
            for j in range(i + 1, len(lines)):
                nextline = lines[j].strip()
                if nextline and not nextline.startswith(";"):
                    parts = nextline.split()
                    if parts:
                        lines[j] = lines[j].replace(parts[0], peg_resname, 1)
                    break
            break
    with open(dst_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"   已复制 PEG 拓扑: {dst_path}")


def _rewrite_itp_moleculetype(itp_path: Path, peg_resname: str) -> None:
    """将 itp 中 [ moleculetype ] 下一行的分子名改为 peg_resname。"""
    with open(itp_path, "r", encoding="utf-8") as f:
        lines = f.read().split("\n")
    for i, line in enumerate(lines):
        if "[ moleculetype ]" in line:
            for j in range(i + 1, len(lines)):
                nextline = lines[j].strip()
                if nextline and not nextline.startswith(";"):
                    parts = nextline.split()
                    if parts and parts[0] != peg_resname:
                        lines[j] = lines[j].replace(parts[0], peg_resname, 1)
                    break
            break
    with open(itp_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def _run_polyply_gen_params(
    peg_length: int,
    head_resname: str,
    output_itp: Path,
    polyply_cfg: dict,
) -> None:
    """用 polyply gen_params 生成 PEG itp（HEAD:1 PEO:N）。参考 version1/pipeline.py ensure_dpeg_itp。"""
    # polyply 头基残基名可单独配置；默认用 config.polyply.head_resname（如 DSPE），
    # 若未配置则回退到 pipe 的 head_resname（如 DPPE）
    head_resname_pol = (polyply_cfg.get("head_resname") or head_resname).strip()
    work_dir = _resolve_path(polyply_cfg.get("work_dir", "default/itps"), VERSION3_DIR)
    head_itp = polyply_cfg.get("head_itp") or f"{head_resname_pol}.itp"
    link_ff = polyply_cfg.get("link_ff", "DSPE_PEO_link.ff")
    lib = polyply_cfg.get("lib", "martini3")
    name = polyply_cfg.get("name", "DPEG")
    if not work_dir.is_dir():
        raise FileNotFoundError(f"polyply work_dir 不存在: {work_dir}")
    head_path = work_dir / head_itp if not Path(head_itp).is_absolute() else Path(head_itp)
    link_path = work_dir / link_ff if not Path(link_ff).is_absolute() else Path(link_ff)
    if not head_path.is_file():
        raise FileNotFoundError(f"polyply 头基 itp 不存在: {head_path}")
    if not link_path.is_file():
        raise FileNotFoundError(f"polyply 连接力场不存在: {link_path}，可从 version1 复制 DSPE_PEO_link.ff 到 {work_dir}")
    cmd = [
        "polyply", "gen_params",
        "-name", name,
        "-seq", f"{head_resname_pol}:1", f"PEO:{peg_length}",
        "-o", str(output_itp.resolve()),
        "-lib", lib,
        "-f", str(head_path.name), str(link_path.name),
    ]
    print(f"   运行: {' '.join(cmd)}")
    subprocess.run(cmd, cwd=str(work_dir), check=True)


def ensure_peg_itp(config: dict, out_dir: Path) -> None:
    """确保 default/itps/{peg_resname}.itp 存在。

    逻辑：
      1) 若 version3/default/itps/{peg_resname}.itp 已存在，直接使用（视为全局库），不再生成；
      2) 否则：
         - 若配置了 dpeg_itp 且文件存在：读入并改 moleculetype 为 peg_resname，写到 default/itps/{peg_resname}.itp；
         - 若 dpeg_itp 不存在或未配置：用 polyply gen_params 在 default/itps 下生成，再改 moleculetype。

    注意：不强制复制到每个 out_dir/itps，而是以 default/itps 为中心统一管理 PEG itp。
    """
    peg_resname = (config.get("peg_resname") or "15PEG").strip()[:5]
    peg_length = int(config.get("peg_length", 15))
    head_resname = (config.get("head_resname") or "DPPE").strip()
    # 全局 itp 库目录：version3/default/itps
    itps_dst = VERSION3_DIR / "default" / "itps"
    itps_dst.mkdir(parents=True, exist_ok=True)
    dst_path = itps_dst / f"{peg_resname}.itp"

    # 1) 已存在则直接返回
    if dst_path.is_file():
        print(f"   检测到已有 PEG itp: {dst_path}，直接使用")
        return

    src = config.get("dpeg_itp")
    if src and str(src).strip():
        src_path = _resolve_path(src.strip(), VERSION3_DIR)
        if src_path.is_file():
            with open(src_path, "r", encoding="utf-8") as f:
                text = f.read()
            lines = text.split("\n")
            for i, line in enumerate(lines):
                if "[ moleculetype ]" in line:
                    for j in range(i + 1, len(lines)):
                        nextline = lines[j].strip()
                        if nextline and not nextline.startswith(";"):
                            parts = nextline.split()
                            if parts:
                                lines[j] = lines[j].replace(parts[0], peg_resname, 1)
                            break
                    break
            with open(dst_path, "w", encoding="utf-8") as f:
                f.write("\n".join(lines))
            print(f"   已复制 PEG 拓扑: {dst_path}")
            return
        print(f"   未找到 dpeg_itp: {src_path}，改用 polyply 生成 …")

    # 无 dpeg_itp 或文件不存在：用 polyply 生成
    polyply_cfg = config.get("polyply") or {}
    try:
        _run_polyply_gen_params(peg_length, head_resname, dst_path, polyply_cfg)
        _rewrite_itp_moleculetype(dst_path, peg_resname)
        print(f"   已用 polyply 生成: {dst_path}")
    except FileNotFoundError as e:
        print(f"   {e}")
        print("   请配置 config.yaml 中 dpeg_itp 指向已有 itp，或确保 polyply 与 default/itps 下 DSPE.itp、DSPE_PEO_link.ff 可用")
    except subprocess.CalledProcessError as e:
        print(f"   polyply gen_params 失败: {e}")
    except Exception as e:
        print(f"   生成 PEG itp 失败: {e}")


def _copy_default_pipe_assets_to_run_dir(
    run_dir: Path,
    output_gro: Path,
    output_top: Path | None,
) -> None:
    """将 `version3/default` 中的资源复制到 run_dir，便于直接运行 default/pipe.sh。

    同时把当前输出文件复制为 default/pipe.sh 所需的文件名：
    - system.gro
    - system.top
    """
    default_dir = VERSION3_DIR / "default"
    out_itps = run_dir / "itps"
    out_mdps = run_dir / "mdps"

    itps_src = default_dir / "itps"
    if itps_src.is_dir():
        shutil.copytree(itps_src, out_itps, dirs_exist_ok=True)

    # 默认脚本里用 mdps/，但目录实际叫 mpds
    mdps_src = default_dir / "mdps"
    if not mdps_src.is_dir():
        mdps_src = default_dir / "mpds"
    if mdps_src.is_dir():
        shutil.copytree(mdps_src, out_mdps, dirs_exist_ok=True)

    pipe_sh_src = default_dir / "pipe.sh"
    if pipe_sh_src.is_file():
        shutil.copy2(pipe_sh_src, run_dir / "pipe.sh")

    # default/pipe.sh 期望 system.gro/system.top/system.ndx
    system_gro = run_dir / "system.gro"
    if output_gro.is_file():
        shutil.copy2(output_gro, system_gro)
    if output_top and output_top.is_file():
        shutil.copy2(output_top, run_dir / "system.top")


def _generate_ndx_for_run_dir(run_dir: Path, config: dict, gro_for_ndx: Path) -> None:
    """使用 version2/ndx.py 生成 system.ndx（真空体系：water/ions 为空）。"""
    ndx_path = run_dir / "system.ndx"

    # 需要包含 PEG 的残基名（最终 gro 里是 peg_resname）
    lipid_resnames = []
    lipids_cfg = config.get("lipids") or {}
    for k in lipids_cfg.keys():
        lipid_resnames.append(str(k))

    peg_resname = (config.get("peg_resname") or "15PEG").strip()[:5]
    lipid_resnames.append(peg_resname)

    gas_resnames = []
    gas = config.get("gas")
    if gas:
        gas_resnames = [str(gas)]

    sys.path.insert(0, str(VERSION2_DIR))
    try:
        from ndx import generate_ndx as do_ndx  # type: ignore
        do_ndx(
            gro_path=str(gro_for_ndx),
            ndx_path=str(ndx_path),
            lipid_resnames=lipid_resnames,
            gas_resnames=gas_resnames or ["O2"],
            water_resnames=[],
            ion_resnames=[],
        )
    finally:
        if str(VERSION2_DIR) in sys.path:
            sys.path.remove(str(VERSION2_DIR))

def _fix_overlaps_in_gro(
    gro_path: Path,
    config: dict,
) -> None:
    """调用 fix_overlap_gro 的逻辑，对 gro 中所有原子做几何松弛，消除 < clash_nm 的过近接触。

    - 默认启用（fix_overlap_enable 默认为 true），可在 config.yaml 设置 fix_overlap_enable: false 关闭；
    - 阈值可通过 fix_overlap_clash_nm / fix_overlap_target_nm / fix_overlap_max_rounds 配置；
    - 直接在原文件上原地更新。
    """
    if not config.get("fix_overlap_enable", True):
        print("3. 跳过 fix_overlap（config.fix_overlap_enable = false）")
        return

    clash_nm = float(config.get("fix_overlap_clash_nm", 0.35))
    target_nm = float(config.get("fix_overlap_target_nm", 0.48))
    max_rounds = int(config.get("fix_overlap_max_rounds", 15))

    sys.path.insert(0, str(VERSION3_DIR))
    try:
        import fix_overlap_gro as fog  # type: ignore
    except ImportError:
        print("3. 未找到 fix_overlap_gro.py，跳过重叠修复")
        sys.path.remove(str(VERSION3_DIR))
        return

    print(f"3. 使用 fix_overlap_gro 对 {gro_path.name} 做几何松弛（clash<{clash_nm:.2f} nm, 目标 {target_nm:.2f} nm）…")
    try:
        title, n_atoms, lines, coords, box = fog.read_gro(str(gro_path))
        print(f"   Atoms: {n_atoms}, box: {box[:3]}")

        # 解析每个原子的 resnr（与 coords/lines 同顺序），避免同一残基内被“推散”导致链断裂/跨 PBC 观感截断
        resids = []
        for l in lines:
            try:
                resids.append(int(l[0:5]))
            except Exception:
                resids.append(-1)

        def _fix_overlaps_cross_residue(coords_nm, box_nm, resids_arr, clash_nm_=0.35, target_nm_=0.48):
            n = len(coords_nm)
            total_moves = 0
            cell_size = 1.0
            nx = max(1, int(box_nm[0] / cell_size))
            ny = max(1, int(box_nm[1] / cell_size))
            nz = max(1, int(box_nm[2] / cell_size))
            cells = {}
            for i in range(n):
                c = coords_nm[i]
                ix = int(c[0] / box_nm[0] * nx) % nx
                iy = int(c[1] / box_nm[1] * ny) % ny
                iz = int(c[2] / box_nm[2] * nz) % nz
                cells.setdefault((ix, iy, iz), []).append(i)

            def pbc_dx(dx, L):
                return dx - round(dx / L) * L

            def pbc_vector(ri, rj, box_):
                return [
                    pbc_dx(rj[0] - ri[0], box_[0]),
                    pbc_dx(rj[1] - ri[1], box_[1]),
                    pbc_dx(rj[2] - ri[2], box_[2]),
                ]

            def norm(v):
                return math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

            for ix in range(nx):
                for iy in range(ny):
                    for iz in range(nz):
                        indices = []
                        for di in (-1, 0, 1):
                            for dj in (-1, 0, 1):
                                for dk in (-1, 0, 1):
                                    key = ((ix + di) % nx, (iy + dj) % ny, (iz + dk) % nz)
                                    if key in cells:
                                        indices.extend(cells[key])
                        for ii in range(len(indices)):
                            for jj in range(ii + 1, len(indices)):
                                i, j = indices[ii], indices[jj]
                                if resids_arr[i] == resids_arr[j] and resids_arr[i] != -1:
                                    continue
                                v = pbc_vector(coords_nm[i], coords_nm[j], box_nm)
                                d = norm(v)
                                if d < 1e-6:
                                    d = 1e-6
                                if d < clash_nm_:
                                    move = (target_nm_ - d) / 2.0
                                    u_ = [v[0]/d, v[1]/d, v[2]/d]
                                    for k in range(3):
                                        coords_nm[i][k] -= move * u_[k]
                                        coords_nm[j][k] += move * u_[k]
                                    for k in range(3):
                                        coords_nm[i][k] = coords_nm[i][k] % box_nm[k]
                                        coords_nm[j][k] = coords_nm[j][k] % box_nm[k]
                                    total_moves += 1
            return total_moves

        for r in range(max_rounds):
            moves = _fix_overlaps_cross_residue(coords, box, resids, clash_nm_=clash_nm, target_nm_=target_nm)
            if moves == 0:
                print(f"   Round {r+1}: no pairs < {clash_nm:.2f} nm, done.")
                break
            print(f"   Round {r+1}: resolved {moves} overlapping pairs")
        tmp_path = gro_path.with_suffix(".gro.tmp_fix")
        fog.write_gro(str(tmp_path), title, n_atoms, lines, coords, box)
        shutil.move(tmp_path, gro_path)
        print(f"   已写回修正后的 {gro_path}")
    except Exception as e:
        print(f"   fix_overlap 过程中出错，保留原始 gro（{gro_path}）：{e}")
    finally:
        if str(VERSION3_DIR) in sys.path:
            sys.path.remove(str(VERSION3_DIR))


def _recenter_resize_and_wrap(
    gro_path: Path,
    config: dict,
) -> None:
    """把 LNB 质心移到盒子中心，并按 PEG 最外半径自动放大盒子。

    目标：
    - 避免体系偏在盒子角落导致 PBC 视觉错乱/链跨边界；
    - 不依赖原始 gro 的 box 尺寸：先以体系自身质心为参考，再决定新 box；

    步骤（单位 nm）：
    1) 以“体系（非溶剂）”的 center of geometry 作为 LNB 质心 c；
    2) 计算 PEG 原子到 c 的最大距离 R_max（若无 PEG 则用全部原子）；
    3) 设置立方盒子边长 L = 2 * R_max + 6 nm；
    4) 将所有坐标平移到盒子中心 (L/2, L/2, L/2)，并 wrap 到 [0, L)；
    5) 写回 gro（由 MDAnalysis 写出，box 同步更新）。
    """
    if not _has_mda:
        print("   跳过自动放大盒子：未安装 MDAnalysis。pip install MDAnalysis")
        return

    u = mda.Universe(str(gro_path))
    peg_resname = (config.get("peg_resname") or "15PEG").strip()[:5]
    solvent_set = {"W", "NA", "CL"}
    non_solvent_ag = u.select_atoms("not resname W NA CL") if u.atoms.n_atoms else u.atoms
    if non_solvent_ag.n_atoms == 0:
        non_solvent_ag = u.atoms

    # 1) LNB 质心：使用非溶剂的 center of geometry
    center_nm = non_solvent_ag.center_of_geometry() / 10.0  # Å -> nm

    # 2) 优先使用 PEG 原子估算 R_max
    ag = u.select_atoms(f"resname {peg_resname}")
    if ag.n_atoms == 0:
        ag = u.atoms
        print(f"   未在 {gro_path.name} 中找到 {peg_resname}，按全部原子估算盒子半径")

    coords_nm = ag.positions / 10.0
    dists = np.linalg.norm(coords_nm - center_nm, axis=1)
    if dists.size == 0:
        print("   自动放大盒子失败：未能获取坐标")
        return
    r_max = float(dists.max())
    L = 2.0 * r_max + 6.0  # nm
    print(f"   PEG 最远距离 R_max = {r_max:.3f} nm，设置立方盒子边长 L = {L:.3f} nm")

    # 3) 平移到盒子中心，并按“分子/残基”为单位做 PBC 处理，避免一条 PEG 被拆成两段
    target_center_nm = np.array([L / 2.0, L / 2.0, L / 2.0], dtype=float)
    shift_nm = target_center_nm - center_nm
    pos_nm = (u.atoms.positions / 10.0) + shift_nm  # nm

    # 先将每个残基内部做最小镜像“成团”（以该残基第一个原子为参考点），再把参考点 wrap 到盒内
    box = np.array([L, L, L], dtype=float)
    for res in u.residues:
        idx = res.atoms.indices
        if idx.size == 0:
            continue
        ref = pos_nm[idx[0]].copy()
        # 将残基内各原子相对 ref 做最小镜像，保证同一残基不跨盒子
        deltas = pos_nm[idx] - ref
        deltas -= np.round(deltas / box) * box
        pos_nm[idx] = ref + deltas
        # 将 ref（以及整段残基）wrap 到 [0, L)
        ref_wrapped = ref % box
        shift_res = ref_wrapped - ref
        pos_nm[idx] += shift_res
        # 进一步确保该残基所有原子都落在 [0, L) 内（整体平移，不破坏相对几何）
        # 用 while 而不是一次 if：处理跨度接近半盒的极端情况，避免出现负坐标导致“被截断”
        for _ in range(3):  # 最多迭代 3 次足够
            mins = pos_nm[idx].min(axis=0)
            maxs = pos_nm[idx].max(axis=0)
            shifted = False
            for ax in range(3):
                if mins[ax] < 0:
                    pos_nm[idx][:, ax] += box[ax]
                    shifted = True
                if maxs[ax] >= box[ax]:
                    pos_nm[idx][:, ax] -= box[ax]
                    shifted = True
            if not shifted:
                break

    u.atoms.positions = pos_nm * 10.0  # Å
    u.dimensions = np.array([L * 10.0, L * 10.0, L * 10.0, 90.0, 90.0, 90.0])

    u.atoms.write(str(gro_path))
    print(f"   已将体系质心移到盒子中心并写回 {gro_path.name}（L={L:.3f} nm）")


def _write_simple_top(
    top_path: Path,
    ordered_residues: list,
    peg_resname: str,
) -> None:
    """根据残基列表写出简易 system_peg.top：自动 include 所需 itp + [ system ] + [ molecules ]。

    - itp 路径统一假定为相对当前 top 的 itps/ 子目录（如 itps/DPPC.itp）；
    - include 顺序参考 version2/build_lnb：martini_v3、ffbond、ions、solu、各脂质、O2、PEG 等；
    - 仅 include 实际存在于 default/itps 下的 itp 文件。
    """
    from collections import Counter
    counts = Counter()
    for r in ordered_residues:
        name = (r.resname or "").strip()
        if name:
            counts[name] += 1
    # 固定顺序优先，其余按名称排序，避免漏写 O2 等导致 grompp 原子数不一致
    order_first = ("DPPC", "DOPC", "DOPS", "CHOL", "DPPE", "DSPE", peg_resname, "O2", "W", "NA", "CL")
    rest = sorted(k for k in counts if k not in order_first)
    # 1) include 段：从默认 itp 库中仅 include 已存在的 itp
    default_itps_dir = VERSION3_DIR / "default" / "itps"
    itp_order = [
        "martini_v3.itp",
        "ffbond.itp",
        "ions.itp",
        "solu.itp",
        "DPPC.itp",
        "DOPS.itp",
        "CHOL.itp",
        "O2.itp",
    ]
    # PEG itp：peg_resname 对应的 itp 文件名（如 45PEG.itp）
    peg_itp_name = f"{peg_resname}.itp"
    itp_order.append(peg_itp_name)
    include_lines = []
    for name in itp_order:
        if (default_itps_dir / name).is_file():
            include_lines.append(f"#include \"itps/{name}\"\n")

    # 2) system 与 molecules 段
    lines = []
    lines.extend(include_lines)
    lines.append("\n; topology after PEG growth, no solvent (vacuum)\n")
    lines.append("[ system ]\n")
    lines.append("LNB PEG\n\n")
    lines.append("[ molecules ]\n")
    lines.append("; name count\n")
    for name in order_first + tuple(rest):
        n = counts.get(name, 0)
        if n > 0:
            lines.append(f"{name:10s} {n:6d}\n")
    top_path.parent.mkdir(parents=True, exist_ok=True)
    with open(top_path, "w", encoding="utf-8") as f:
        f.writelines(lines)
    print(f"   已写 {top_path}")


def load_config(config_path: Path | None = None) -> dict:
    if config_path is None:
        config_path = VERSION3_DIR / "config.yaml"
    if not config_path.is_file():
        return {}
    with open(config_path, "r", encoding="utf-8") as f:
        raw = f.read()
    if yaml:
        return yaml.safe_load(raw) or {}
    # 简单解析
    data = {}
    for line in raw.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or ":" not in line:
            continue
        k, _, v = line.partition(":")
        k, v = k.strip(), v.strip()
        if v.lower() == "true":
            data[k] = True
        elif v.lower() == "false":
            data[k] = False
        elif v.isdigit():
            data[k] = int(v)
        else:
            try:
                data[k] = float(v)
            except ValueError:
                data[k] = v
    return data


def main():
    config_path = VERSION3_DIR / "config.yaml"
    if len(sys.argv) > 1:
        config_path = Path(sys.argv[1])
    config = load_config(config_path)
    if not config:
        print("未找到或无法解析 config.yaml")
        sys.exit(1)

    base = VERSION3_DIR
    # 所有生成物统一输出到 run_dir
    run_dir = _resolve_path(config.get("run_dir") or config.get("output_dir") or "output", VERSION3_DIR)
    run_dir.mkdir(parents=True, exist_ok=True)

    input_gro = base / config.get("input_gro", "system.gro")

    out_gro_name = config.get("output_gro", "system_peg.gro")
    output_gro = run_dir / out_gro_name if not Path(str(out_gro_name)).is_absolute() else Path(str(out_gro_name))

    out_top = config.get("output_top")
    if out_top:
        output_top = run_dir / out_top if not Path(str(out_top)).is_absolute() else Path(str(out_top))
    else:
        output_top = None

    if not input_gro.is_file():
        print(f"输入 GRO 不存在: {input_gro}")
        sys.exit(1)

    print("1. 在 DPPE 的 NH3 上沿球心法线向外生长 PEG …")
    gro_with_peg = output_gro.parent / (output_gro.stem + "_with_peg.gro")
    if gro_with_peg == output_gro:
        gro_with_peg = output_gro.parent / "system_with_peg.gro"
    grow_peg_on_dppe(input_gro, gro_with_peg, config)

    print("2. 移除全部水与溶剂离子（用于真空模拟）…")
    remove_all_solvent(gro_with_peg, config, output_gro, output_top)
    if gro_with_peg != output_gro and gro_with_peg.exists():
        gro_with_peg.unlink()

    # 3. 自动做一次几何松弛（fix_overlap），缓解 PEG 链与其他 beads 的过近接触
    _fix_overlaps_in_gro(output_gro, config)

    # 4. 质心归中 + 按 PEG 最外半径自动放大盒子：L = 2 * R_max + 6 nm，并 wrap 到盒内
    print("4. 质心归中并按 PEG 最外半径自动调整盒子大小 …")
    _recenter_resize_and_wrap(output_gro, config)

    # 5. 确保 PEG itp 存在（有 dpeg_itp 则复用，否则用 polyply 生成到 default/itps）
    out_dir = output_gro.parent
    print("5. 确保 PEG itp 存在（有 dpeg_itp 则复制，否则 polyply 生成）…")
    ensure_peg_itp(config, out_dir)

    print("6. 复制默认模拟资源并生成 system.ndx …")
    _copy_default_pipe_assets_to_run_dir(out_dir, output_gro, output_top)
    _generate_ndx_for_run_dir(out_dir, config, out_dir / "system.gro")

    print("完成。")


if __name__ == "__main__":
    main()
