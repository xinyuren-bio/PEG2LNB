#!/usr/bin/env python3
"""
根据 YAML 配置：对“干体系”执行
1) gmx solvate 填充 Martini 水（W）
2) MDAnalysis 删除气泡质心半径 bubble_radius_nm 内的指定残基（建议仅 W；此时尚无盐离子）
3) gmx genion -neutral：用水分子置换为 NA/CL，使体系净电荷为 0
4) gmx genion -conc：在中和后再按目标盐浓度补加等量 NA/CL
5) 按最终 gro 回写 top 中 W/NA/CL 的 [ molecules ] 计数；
   genion 先 -neutral 再 -conc 时，.gro 离子区常为「残基名 NA（成盐 Na⁺）→ CL → ION/NA（中和 Na⁺）」，
   与 itp 中分子类型名均为 NA，但 top 须写两行 NA（分段计数），顺序为：
   W → NA(成盐) → CL → NA(中和)。

输入期望：
- in_gro 与 in_top 对应，且 in_top 至少能让 gmx solvate/genion 正常写入/修改 molecules
输出：
- out_solvated_gro/out_solvated_top：仅溶剂化（检查点）
- out_removed_gro/out_removed_top：删泡内溶剂并中和后的最终体系
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any


def _load_yaml(path: Path) -> dict[str, Any]:
    try:
        import yaml  # type: ignore
    except Exception as e:  # pragma: no cover
        raise RuntimeError(
            "需要 PyYAML：pip install pyyaml 才能解析 yaml 配置文件。"
        ) from e

    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def _run(cmd: list[str], *, cwd: Path | None = None, stdin_text: str | None = None) -> None:
    print("\n[CMD] " + " ".join(cmd))
    if stdin_text is not None:
        print("[STDIN]\n" + stdin_text.rstrip("\n"))
    p = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        input=(stdin_text.encode("utf-8") if stdin_text is not None else None),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    out = p.stdout.decode("utf-8", errors="replace") if p.stdout else ""
    if out.strip():
        print("[gmx output]\n" + "\n".join(out.splitlines()[-80:]))
    if p.returncode != 0:
        raise RuntimeError(f"命令执行失败，returncode={p.returncode}: {' '.join(cmd)}")


def _ensure_include_lines(top_text: str, include_lines: list[str]) -> str:
    """
    如果 top_text 中没有指定 include 行（子串匹配），则插入到所有 #include 之后。
    """
    missing = [inc for inc in include_lines if inc not in top_text]
    if not missing:
        return top_text

    lines = top_text.splitlines(keepends=True)
    # 找到最后一条 include 的位置（只在文件头部范围内，避免把 molecules 里注释也算进去）
    last_inc_idx = -1
    for i, line in enumerate(lines):
        if i > 200:
            break
        if line.strip().startswith("#include"):
            last_inc_idx = i
    insert_at = last_inc_idx + 1 if last_inc_idx >= 0 else 0
    for inc in missing:
        lines.insert(insert_at, f"{inc}\n")
        insert_at += 1
    return "".join(lines)


def _find_section_bounds(text_lines: list[str], section_name: str) -> tuple[int, int] | None:
    """
    返回 [start,end) 范围：start 行为 [ section_name ] 那行，end 为下一个 [ ... ] 起始或 EOF。
    """
    sec_re = re.compile(rf"^\s*\[\s*{re.escape(section_name)}\s*\]\s*$", re.IGNORECASE)
    start = None
    for i, line in enumerate(text_lines):
        if sec_re.match(line):
            start = i
            break
    if start is None:
        return None
    end = len(text_lines)
    for j in range(start + 1, len(text_lines)):
        if text_lines[j].lstrip().startswith("[") and re.match(r"^\s*\[.*\]\s*$", text_lines[j]):
            end = j
            break
    return start, end


def _update_molecules_counts_in_top(top_path: Path, counts: dict[str, int]) -> None:
    """
    更新 [ molecules ] 段中指定分子名的计数；若缺失则插入。
    counts: {"W":123,...}
    """
    text = top_path.read_text(encoding="utf-8")
    lines = text.splitlines(keepends=True)
    bounds = _find_section_bounds(lines, "molecules")
    if not bounds:
        raise RuntimeError(f"未找到 top 中 [ molecules ] 段: {top_path}")
    start, end = bounds

    # 记录每个分子是否已出现并被替换
    updated = {k.upper(): False for k in counts.keys()}
    wanted = {k.upper() for k in counts.keys()}

    new_lines: list[str] = []
    new_lines.extend(lines[:start + 1])

    # 保留 [ molecules ] 行之后的行，直到 end，按需替换
    for i in range(start + 1, end):
        line = lines[i]
        stripped = line.strip()
        if not stripped or stripped.startswith(";") or stripped.startswith("#"):
            new_lines.append(line)
            continue
        parts = stripped.split()
        if not parts:
            new_lines.append(line)
            continue
        mol = parts[0].strip().upper()
        if mol in wanted:
            # genion -neutral 有时会在 [ molecules ] 里追加第二行同名离子；
            # 只保留第一行并写入与 gro 一致的计数，避免两行同名导致 GROMACS 误算总数。
            if updated.get(mol):
                continue
            cnt = counts.get(mol, 0)
            # 兼容多空格；保持名字列宽不强制
            new_lines.append(f"{parts[0]:<16}{cnt:>6d}\n")
            updated[mol] = True
        else:
            new_lines.append(line)

    # 插入缺失分子到段尾（end 前），放在最后非注释行之后即可
    insert_at = len(new_lines)
    for mol, ok in updated.items():
        if not ok:
            insert_at = len(new_lines)
            # 用标准格式插入
            new_lines.insert(insert_at, f"{mol:<16}{counts[mol]:>6d}\n")

    top_path.write_text("".join(new_lines), encoding="utf-8")


def _update_molecules_w_na_cl_na_order(
    top_path: Path,
    *,
    w: int,
    na_salt: int,
    cl: int,
    na_neutral: int,
) -> None:
    """
    与 genion 先 -neutral 再 -conc 的典型 .gro 一致：成盐段残基名常为 NA，中和段残基名为 ION、
    原子名为 NA；两段在拓扑里均为分子类型 NA，须分两行，且顺序为
    W → NA(成盐) → CL → NA(中和)。
    """
    text = top_path.read_text(encoding="utf-8")
    lines = text.splitlines(keepends=True)
    bounds = _find_section_bounds(lines, "molecules")
    if not bounds:
        raise RuntimeError(f"未找到 top 中 [ molecules ] 段: {top_path}")
    start, end = bounds

    new_lines: list[str] = []
    new_lines.extend(lines[: start + 1])

    inserted = False
    for i in range(start + 1, end):
        line = lines[i]
        stripped = line.strip()
        if not stripped or stripped.startswith(";") or stripped.startswith("#"):
            new_lines.append(line)
            continue
        parts = stripped.split()
        mol = parts[0].strip().upper()
        if mol in {"W", "NA", "CL"}:
            if not inserted:
                new_lines.append(f"{'W':<16}{w:>6d}\n")
                new_lines.append(f"{'NA':<16}{na_salt:>6d}\n")
                new_lines.append(f"{'CL':<16}{cl:>6d}\n")
                new_lines.append(f"{'NA':<16}{na_neutral:>6d}\n")
                inserted = True
            continue
        new_lines.append(line)

    if not inserted:
        raise RuntimeError(
            f"[ molecules ] 中未找到 W/NA/CL 行，无法写入离子计数: {top_path}"
        )

    new_lines.extend(lines[end:])
    top_path.write_text("".join(new_lines), encoding="utf-8")


def _extract_molecule_names_from_top(top_text: str) -> list[str]:
    """
    提取 [ molecules ] 段中声明过的分子名（按出现顺序去重）。
    """
    lines = top_text.splitlines(keepends=False)
    bounds = _find_section_bounds([l + "\n" for l in lines], "molecules")
    if not bounds:
        return []
    start, end = bounds
    names: list[str] = []
    seen: set[str] = set()
    for i in range(start + 1, end):
        stripped = lines[i].strip()
        if not stripped or stripped.startswith(";") or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if not parts:
            continue
        name = parts[0].strip()
        key = name.upper()
        if key not in seen:
            seen.add(key)
            names.append(name)
    return names


def _count_residues_by_resname(u: Any, resnames: set[str]) -> dict[str, int]:
    resnames_upper = {r.upper() for r in resnames}
    out: dict[str, int] = {r.upper(): 0 for r in resnames_upper}
    for r in u.residues:
        rn = (r.resname or "").strip().upper()
        if rn in out:
            out[rn] += 1
    return out


def _count_solvent_and_ions_for_top(u: Any) -> dict[str, int]:
    """
    统计用于回写 [ molecules ] 的 W/NA/CL 计数。

    Martini `ions.itp` 里的 NA/CL 分子常写成：
    - moleculetype: NA / CL
    - residu(resname): ION
    - atomname: NA / CL
    因此这里对离子按残基逐个识别：
    - W: resname == W
    - NA/CL: 优先识别 resname == NA/CL；否则兼容 resname == ION 且首原子名为 NA/CL
    """
    counts = {"W": 0, "NA": 0, "CL": 0}
    for r in u.residues:
        rn = (r.resname or "").strip().upper()
        if rn == "W":
            counts["W"] += 1
            continue
        if rn in {"NA", "CL"}:
            counts[rn] += 1
            continue
        if rn == "ION" and len(r.atoms) > 0:
            atom_name = (r.atoms[0].name or "").strip().upper()
            if atom_name in {"NA", "CL"}:
                counts[atom_name] += 1
    return counts


def _min_image_dist_nm(p_nm: Any, center_nm: Any, box_nm: Any) -> float:
    """
    p_nm/center_nm/box_nm 都是 (3,) 的 nm 标量数组
    """
    import numpy as np

    d = np.array(p_nm, dtype=float) - np.array(center_nm, dtype=float)
    d = d - np.round(d / np.array(box_nm, dtype=float)) * np.array(box_nm, dtype=float)
    return float(np.linalg.norm(d))


def _write_gro_manual(
    residues_kept: list[Any],
    out_gro: Path,
    box_nm: tuple[float, float, float],
    title: str,
    *,
    solvent_atom_serial_start: int = 1,
):
    """
    residues_kept: 按原始顺序排列的 Residue 列表
    用固定列宽写 gro，保证 grompp 可读。
    """
    fmt = "{resnr:5d}{resname:<5}{atomname:<5}{atomnr:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n"

    out_gro.parent.mkdir(parents=True, exist_ok=True)

    # 重新编号 resid/residue、全局 atomnr
    atomnr = solvent_atom_serial_start
    resnr = 0
    total_atoms = sum(len(r.atoms) for r in residues_kept)

    lines: list[str] = []
    lines.append(title.rstrip() + "\n")
    lines.append(f"{total_atoms:5d}\n")

    for r in residues_kept:
        resnr += 1
        # GRO 的 resnr/atomnr 都是 5 位，超过需要 wrap
        resnr_out = ((resnr - 1) % 99999) + 1
        rname = (r.resname or "").strip()[:5]
        for a in r.atoms:
            atomname = (a.name or "").strip()[:5]
            x, y, z = a.position  # Å
            atomnr_out = ((atomnr - 1) % 99999) + 1
            lines.append(
                fmt.format(
                    resnr=resnr_out,
                    resname=rname,
                    atomname=atomname,
                    atomnr=atomnr_out,
                    x=x / 10.0,
                    y=y / 10.0,
                    z=z / 10.0,
                )
            )
            atomnr += 1

    lines.append(f"{box_nm[0]:10.5f}{box_nm[1]:10.5f}{box_nm[2]:10.5f}\n")
    out_gro.write_text("".join(lines), encoding="utf-8")


def _parse_group_input(group: str) -> str:
    # genion prompt 支持输入 group 名/序号（不同版本不同）；这里不做强制转换
    return group.strip()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=str, required=True, help="YAML 配置路径")
    args = ap.parse_args()

    cfg_path = Path(args.config).expanduser().resolve()
    if not cfg_path.is_file():
        raise FileNotFoundError(f"配置文件不存在: {cfg_path}")

    cfg = _load_yaml(cfg_path)

    repo_root = Path(__file__).resolve().parent

    def _resolve_repo(p: str) -> Path:
        path = Path(str(p).strip().strip('"').strip("'"))
        if not path.is_absolute():
            path = repo_root / path
        return path.resolve()

    # --- 输入输出 ---
    in_gro = _resolve_repo(cfg["in_gro"])
    in_top = _resolve_repo(cfg["in_top"])
    output_dir = _resolve_repo(cfg.get("output_dir", str(in_gro.parent)))
    output_dir.mkdir(parents=True, exist_ok=True)

    solvent_cs_gro = Path(
        cfg.get(
            "solvent_cs_gro",
            "default/shared_assets/martini_W.gro",
        )
    )
    if not solvent_cs_gro.is_absolute():
        solvent_cs_gro = repo_root / solvent_cs_gro
    if not solvent_cs_gro.is_file():
        raise FileNotFoundError(f"solvent_cs_gro 不存在: {solvent_cs_gro}")

    bubble_radius_nm = float(cfg.get("bubble_radius_nm", 12.0))
    salt_conc_m = float(cfg.get("salt_conc_m", 0.15))

    # 可选：在溶剂化之前重设盒子尺寸，并把“体系中心”平移到新盒子中心。
    # 用于修复“修改了 gro 的 box 但坐标未平移”的情况。
    target_box_nm = cfg.get("target_box_nm", None)
    in_gro_for_solvation = in_gro
    if target_box_nm is not None:
        try:
            import MDAnalysis as mda  # type: ignore
            import numpy as np  # type: ignore
        except Exception as e:  # pragma: no cover
            raise RuntimeError("需要 MDAnalysis：pip install MDAnalysis") from e

        if isinstance(target_box_nm, (int, float)):
            target_box_tuple = (float(target_box_nm), float(target_box_nm), float(target_box_nm))
        elif isinstance(target_box_nm, (list, tuple)) and len(target_box_nm) >= 3:
            target_box_tuple = (float(target_box_nm[0]), float(target_box_nm[1]), float(target_box_nm[2]))
        else:
            raise ValueError(f"target_box_nm 的值不合法: {target_box_nm}")

        u_in = mda.Universe(str(in_gro))
        if u_in.dimensions is None or len(u_in.dimensions) < 3:
            raise RuntimeError("MDAnalysis 读取输入 gro 的 box 失败（u_in.dimensions 为空）。")

        # 体系中心：用输入 gro 中所有原子的 center_of_geometry
        # （干体系里不含 W/NA/CL，因此天然对应“体系中心”）
        old_center_nm = u_in.atoms.center_of_geometry() / 10.0  # Å -> nm
        new_center_nm = np.array(target_box_tuple, dtype=float) / 2.0
        shift_nm = new_center_nm - old_center_nm

        pos_nm = u_in.atoms.positions / 10.0 + shift_nm
        L = np.array(target_box_tuple, dtype=float)
        pos_nm = np.mod(pos_nm, L)  # wrap 到新盒子 [0,L)
        u_in.atoms.positions = pos_nm * 10.0
        u_in.dimensions = np.array(
            [target_box_tuple[0] * 10.0, target_box_tuple[1] * 10.0, target_box_tuple[2] * 10.0, 90.0, 90.0, 90.0]
        )

        in_gro_for_solvation = output_dir / "in_gro_reboxed.gro"
        u_in.atoms.write(str(in_gro_for_solvation))

        print("\n[Rebox before solvate]")
        print(f"Target box (nm) = {target_box_tuple}")
        print(f"Old system center (nm) = {old_center_nm.tolist()}")
        print(f"New box center (nm) = {new_center_nm.tolist()}")
        print(f"Shift applied (nm) = {shift_nm.tolist()}")
        print(f"Reboxed gro written to: {in_gro_for_solvation}")

    # 质心选择：只用 DPPC/DOPS/CHOL（不含 PEG，不含 O2）
    center_resnames = cfg.get("center_resnames", ["DPPC", "DOPS", "CHOL"])
    center_resnames = [str(x).upper() for x in center_resnames]

    # 删除对象：默认只删泡内 W（加盐前无 NA/CL；若列表含 NA/CL 仅当 gro 里已存在时生效）
    remove_resnames = cfg.get("remove_resnames_inside_radius", ["W"])
    remove_resnames = [str(x).upper() for x in remove_resnames]
    remove_set = set(remove_resnames)

    # genion 替换组：交互提示处输入该字符串（名字或编号）
    # gmx genion 交互时通常会要求选择“要替换的溶剂组”，在 Martini 水体系里组名一般叫 `W`
    genion_solvent_group = str(cfg.get("genion_solvent_group", "W"))
    genion_solvent_group = _parse_group_input(genion_solvent_group)
    gmx_exec = str(cfg.get("gmx_exec", "gmx"))

    # genion 的 mdp（只用于 grompp 生成 tpr）
    genion_mdp = cfg.get("genion_mdp")
    if genion_mdp:
        genion_mdp_path = Path(genion_mdp).expanduser().resolve()
    else:
        genion_mdp_path = (repo_root / "build_output" / "genion.mdp").resolve()
    if not genion_mdp_path.is_file():
        raise FileNotFoundError(f"genion_mdp 不存在: {genion_mdp_path}")

    # 输出文件名
    out_solvated_gro = output_dir / cfg.get("out_solvated_gro", "out_solvated.gro")
    out_solvated_top = output_dir / cfg.get("out_solvated_top", "out_solvated.top")
    out_removed_gro = output_dir / cfg.get("out_removed_gro", "out_removed.gro")
    out_removed_top = output_dir / cfg.get("out_removed_top", "out_removed.top")

    # --- 复制输入 top -> temp top（让 solvate/genion 在临时文件上工作） ---
    tmp_top = output_dir / "tmp_topol.top"
    # 你的 top 通常用相对 include：#include "itps/xxx.itp"
    # gmx 在不同工作目录下对 include 的搜索行为可能不同，因此直接把 in_top 同级的 itps 拷到 output_dir/itps。
    src_itps_dir = in_top.parent / "itps"
    dst_itps_dir = output_dir / "itps"
    if src_itps_dir.is_dir():
        shutil.copytree(src_itps_dir, dst_itps_dir, dirs_exist_ok=True)
    shutil.copy2(in_top, tmp_top)
    tmp_top_text = tmp_top.read_text(encoding="utf-8")
    # 避免 dry top 缺 include 导致 ions/solvent 类型未定义
    tmp_top_text = _ensure_include_lines(
        tmp_top_text,
        [
            '#include "itps/solu.itp"',
            '#include "itps/ions.itp"',
        ],
    )
    # 自动补齐 [ molecules ] 中已出现、且 itps 目录里存在的分子 include（例如 18PEG.itp）
    auto_include_lines: list[str] = []
    for mol_name in _extract_molecule_names_from_top(tmp_top_text):
        itp_path = dst_itps_dir / f"{mol_name}.itp"
        if itp_path.is_file():
            auto_include_lines.append(f'#include "itps/{mol_name}.itp"')
    if auto_include_lines:
        tmp_top_text = _ensure_include_lines(tmp_top_text, auto_include_lines)
    tmp_top.write_text(tmp_top_text, encoding="utf-8")

    # --- Step 1: solvate ---
    # gmx solvate 会更新 -p 指定的 topol.top（molecules 计数）
    _run(
        [
            gmx_exec,
            "solvate",
            "-cp",
            str(in_gro_for_solvation),
            "-cs",
            str(solvent_cs_gro),
            "-o",
            str(out_solvated_gro),
            "-p",
            str(tmp_top),
        ],
        cwd=output_dir,
    )
    # solvate 已更新 tmp_top；此处把它复制成 out_solvated_top 便于检查
    shutil.copy2(tmp_top, out_solvated_top)

    # --- Step 2: MDAnalysis 删除气泡内指定残基（建议仅 W）---
    try:
        import MDAnalysis as mda  # type: ignore
        import numpy as np  # type: ignore
    except Exception as e:  # pragma: no cover
        raise RuntimeError("需要安装 MDAnalysis：pip install MDAnalysis") from e

    u = mda.Universe(str(out_solvated_gro))
    dims = u.dimensions
    if dims is None or len(dims) < 3:
        raise RuntimeError("MDAnalysis 读取 gro 的 box 失败（u.dimensions 为空）。")
    box_nm = (float(dims[0]) / 10.0, float(dims[1]) / 10.0, float(dims[2]) / 10.0)

    # 质心 center：只用指定脂质残基名的所有原子（不含 PEG，不含 O2）
    center_sel = " or ".join(f"resname {r}" for r in center_resnames)
    center_atoms = u.select_atoms(center_sel)
    if center_atoms.n_atoms == 0:
        raise RuntimeError(f"未选中质心原子：center_resnames={center_resnames}，选择表达式={center_sel}")
    center_nm = np.array(center_atoms.centroid()) / 10.0  # Å -> nm（centroid 与 COM 都行，这里用 centroid 更稳）

    # dist<=bubble_radius_nm 内的 remove_set 残基删除
    target_residues = [
        r for r in u.residues
        if (r.resname or "").strip().upper() in remove_set
    ]

    # 对于 W（单珠水/一个原子），r.atoms[0] 即坐标；离子同理
    kept_residues: list[Any] = []
    removed = 0
    for r in u.residues:
        rn = (r.resname or "").strip().upper()
        if rn not in remove_set:
            kept_residues.append(r)
            continue

        # 取 residue 第一个原子作为位置代表（对单珠/离子足够；若有多原子水模型也能适配：此处统一用原子0）
        # 你的 martini 水坐标文件（martini_W.gro）是单珠水，符合此假设。
        if len(r.atoms) == 0:
            kept_residues.append(r)
            continue
        pos_atom = r.atoms[0].position / 10.0  # Å -> nm
        d = _min_image_dist_nm(pos_atom, center_nm, box_nm)
        if d <= bubble_radius_nm:
            removed += 1
        else:
            kept_residues.append(r)

    # 输出 gro
    title = getattr(u.trajectory, "title", None) or "out_removed"
    _write_gro_manual(
        residues_kept=kept_residues,
        out_gro=out_removed_gro,
        box_nm=box_nm,
        title=str(title),
    )

    # --- 统计删掉多少，并写出 updated top ---
    # 直接基于内存 kept_residues 统计，避免再次从 out_removed_gro 读入触发 GRO 解析问题
    counts_kept_upper: dict[str, int] = {k: 0 for k in remove_set}
    for r in kept_residues:
        rn = (r.resname or "").strip().upper()
        if rn in counts_kept_upper:
            counts_kept_upper[rn] += 1

    # 删泡后 top：与 solvate 后的 tmp_top 对齐 W（及 remove_set 中其它项计数）
    shutil.copy2(tmp_top, out_removed_top)
    _update_molecules_counts_in_top(out_removed_top, counts_kept_upper)

    before_counts = _count_residues_by_resname(u, remove_set)
    after_bubble_counts = counts_kept_upper

    def _fmt_counts(d: dict[str, int]) -> str:
        return ", ".join([f"{k}={d.get(k,0)}" for k in sorted(d.keys())])

    print("\n==== Summary (after bubble removal) ====")
    print(f"Bubble center (nm) = {center_nm.tolist()}")
    print(f"Bubble radius (nm) = {bubble_radius_nm}")
    print(f"Removed residues inside bubble (remove_set): {removed}")
    print(f"Before bubble removal: {_fmt_counts({k.upper(): int(v) for k,v in before_counts.items()})}")
    print(f"After  bubble removal: {_fmt_counts({k.upper(): int(v) for k,v in after_bubble_counts.items()})}")
    print(f"out_solvated_gro: {out_solvated_gro}")
    print(f"out_solvated_top: {out_solvated_top}")

    # --- Step 3: genion -neutral（置换 W 为 NA/CL，净电荷 -> 0）---
    ions_tpr = output_dir / "ions.tpr"
    _run(
        [
            gmx_exec,
            "grompp",
            "-f",
            str(genion_mdp_path),
            "-c",
            str(out_removed_gro),
            "-p",
            str(out_removed_top),
            "-o",
            str(ions_tpr),
            "-maxwarn",
            "1",
        ],
        cwd=output_dir,
    )
    _run(
        [
            gmx_exec,
            "genion",
            "-s",
            str(ions_tpr),
            "-o",
            str(out_removed_gro),
            "-p",
            str(out_removed_top),
            "-pname",
            "NA",
            "-nname",
            "CL",
            "-neutral",
        ],
        cwd=output_dir,
        stdin_text=f"{genion_solvent_group}\n{genion_solvent_group}\n",
    )
    u_final = mda.Universe(str(out_removed_gro))
    final_counts = _count_solvent_and_ions_for_top(u_final)
    _update_molecules_counts_in_top(
        out_removed_top, {k.upper(): int(v) for k, v in final_counts.items()}
    )
    na_after_neutral = int(final_counts["NA"])

    print("\n==== Summary (after genion -neutral) ====")
    print(
        "W/NA/CL: "
        + _fmt_counts({k.upper(): int(v) for k, v in final_counts.items()})
    )

    # --- Step 4: genion -conc（中和后再补到目标盐浓度）---
    if salt_conc_m > 0:
        _run(
            [
                gmx_exec,
                "grompp",
                "-f",
                str(genion_mdp_path),
                "-c",
                str(out_removed_gro),
                "-p",
                str(out_removed_top),
                "-o",
                str(ions_tpr),
                "-maxwarn",
                "1",
            ],
            cwd=output_dir,
        )
        _run(
            [
                gmx_exec,
                "genion",
                "-s",
                str(ions_tpr),
                "-o",
                str(out_removed_gro),
                "-p",
                str(out_removed_top),
                "-pname",
                "NA",
                "-nname",
                "CL",
                "-conc",
                str(salt_conc_m),
            ],
            cwd=output_dir,
            stdin_text=f"{genion_solvent_group}\n{genion_solvent_group}\n",
        )
        u_final = mda.Universe(str(out_removed_gro))
        final_counts = _count_solvent_and_ions_for_top(u_final)
        na_salt = int(final_counts["NA"]) - na_after_neutral
        _update_molecules_w_na_cl_na_order(
            out_removed_top,
            w=int(final_counts["W"]),
            na_salt=na_salt,
            cl=int(final_counts["CL"]),
            na_neutral=na_after_neutral,
        )
        print(f"\n==== Summary (after genion -conc {salt_conc_m} M) ====")
        print(
            "W/NA/CL: "
            + _fmt_counts({k.upper(): int(v) for k, v in final_counts.items()})
        )

    print(f"out_removed_gro: {out_removed_gro}")
    print(f"out_removed_top: {out_removed_top}")

    print("\n==== Validation suggestions ====")
    print("1) 用你的 mdp 跑一次 grompp，确认 gro/top 原子数与 molecules 计数一致。")
    print("   例：")
    print(f"     {gmx_exec} grompp -f <your_mdp.mdp> -c {out_removed_gro} -p {out_removed_top} -maxwarn 1")
    print("2) 重点检查 grompp 输出中的：")
    print("   - topology/molecule count mismatch 类错误（这通常表示 top 没跟上 gro 的删改）")
    print("   - Total charge 是否接近 0（中和后再加等量盐，理论上仍应接近 0）")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:  # pragma: no cover
        print("\nInterrupted.")
        raise

