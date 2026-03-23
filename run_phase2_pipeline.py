#!/usr/bin/env python
# Phase2 一键运行：改盒子 -> solvate+genion -> 删气泡内 W/NA/CL -> 生成 ndx/pipe.sh -> 运行含水含盐模拟

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path


def _load_yaml_min(path: Path) -> dict:
    try:
        import yaml  # type: ignore

        return yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except Exception:
        # 简易解析（仅用于 demo 配置字段）
        data: dict = {}
        stack: list[tuple[int, dict]] = [(-1, data)]
        for raw in path.read_text(encoding="utf-8").splitlines():
            line = raw.rstrip()
            if not line or line.strip().startswith("#"):
                continue
            indent = len(line) - len(line.lstrip(" "))
            if ":" not in line:
                continue
            key, _, val = line.strip().partition(":")
            key = key.strip()
            val = val.strip()
            while stack and indent <= stack[-1][0]:
                stack.pop()
            cur = stack[-1][1] if stack else data
            if val == "":
                cur[key] = {}
                stack.append((indent, cur[key]))
            else:
                if val.lower() in ("true", "false"):
                    cur[key] = val.lower() == "true"
                else:
                    try:
                        cur[key] = float(val)
                    except Exception:
                        try:
                            cur[key] = int(val)
                        except Exception:
                            cur[key] = val.strip('"').strip("'")
        return data


def _parse_resnames_from_gro(gro_path: Path) -> set[str]:
    resnames: set[str] = set()
    lines = gro_path.read_text(encoding="utf-8").splitlines()
    n_atoms = int(lines[1].strip())
    for i in range(2, 2 + n_atoms):
        line = lines[i]
        resname = line[5:10].strip()
        if resname:
            resnames.add(resname.upper())
    return resnames


def main():
    repo_root = Path(__file__).resolve().parent  # PEG2LNB
    cfg_path = repo_root / "configs" / "phase2.yaml"
    cfg = _load_yaml_min(cfg_path)

    phase1_run_dir = repo_root / cfg["phase1_run_dir"]
    phase2_run_dir = repo_root / cfg["phase2_run_dir"]
    phase2_run_dir.mkdir(parents=True, exist_ok=True)

    dry_gro = None
    for cand in cfg.get("dry_gro_candidates", []):
        p = phase1_run_dir / cand
        if p.is_file():
            dry_gro = p
            break
    if dry_gro is None:
        # 退一步：直接用 system.gro
        dry_gro = phase1_run_dir / "system.gro"
    dry_top = phase1_run_dir / "system.top"

    # 临时 phase2 配置给 add_solvate_ions_remove_bubble.py 使用
    tmp_cfg = phase2_run_dir / "phase2_add_solvate_cfg.yaml"
    tmp_cfg_text = f"""in_gro: "{dry_gro}"
in_top: "{dry_top}"
output_dir: "{phase2_run_dir}"

salt_conc_m: {cfg.get('salt_conc_m', 0.15)}
bubble_radius_nm: {cfg.get('bubble_radius_nm', 12.0)}

center_resnames: {cfg.get('center_resnames', ['DPPC','DOPS','CHOL'])}
remove_resnames_inside_radius: {cfg.get('remove_resnames_inside_radius', ['W','NA','CL'])}

target_box_nm: {cfg.get('target_box_nm', 30.0)}

solvent_cs_gro: "{cfg.get('solvent_cs_gro', 'version1/shared_assets/martini_W.gro')}"
genion_solvent_group: "{cfg.get('genion_solvent_group', 'W')}"
gas: "{cfg.get('gas','O2')}"

genion_mdp: "{repo_root / 'build_output' / 'genion.mdp'}"

out_removed_gro: "system.gro"
out_removed_top: "system.top"
"""
    tmp_cfg.write_text(tmp_cfg_text, encoding="utf-8")

    print("\n[Phase2] Solvate + add ions + remove bubble solvent/ions ...")
    subprocess.run(
        [sys.executable, str(repo_root / "add_solvate_ions_remove_bubble.py"), "--config", str(tmp_cfg)],
        check=True,
        cwd=str(repo_root),
    )

    # 复制默认 mdps/itps/pipe.sh 到 phase2_run_dir（用于后续模拟）
    default_dir = repo_root / "default"
    # itps
    shutil.copytree(default_dir / "itps", phase2_run_dir / "itps", dirs_exist_ok=True)
    # mdps
    mpds_src = default_dir / "mpds"
    mdps_src = default_dir / "mdps"
    if mdps_src.is_dir():
        shutil.copytree(mdps_src, phase2_run_dir / "mdps", dirs_exist_ok=True)
    elif mpds_src.is_dir():
        shutil.copytree(mpds_src, phase2_run_dir / "mdps", dirs_exist_ok=True)
    # pipe.sh
    shutil.copy2(default_dir / "pipe.sh", phase2_run_dir / "pipe.sh")

    # 生成 system.ndx（Phase2 需要 tc-grps: lnb_layer/lnb_gas/solute）
    print("\n[Phase2] Generating system.ndx ...")
    gro_for_ndx = phase2_run_dir / "system.gro"
    ndx_path = phase2_run_dir / "system.ndx"
    gas_resname = str(cfg.get("gas", "O2")).upper()
    water_resnames = ["W"]
    ion_resnames = ["NA", "CL"]

    all_resnames = _parse_resnames_from_gro(gro_for_ndx)
    lipid_resnames = sorted([r for r in all_resnames if r not in set([gas_resname] + water_resnames + ion_resnames)])

    sys.path.insert(0, str(repo_root))
    try:
        from ndx import generate_ndx  # type: ignore

        generate_ndx(
            gro_path=str(gro_for_ndx),
            ndx_path=str(ndx_path),
            lipid_resnames=lipid_resnames,
            gas_resnames=[gas_resname],
            water_resnames=water_resnames,
            ion_resnames=ion_resnames,
        )
    finally:
        sys.path.remove(str(repo_root))

    # 不自动跑 MD：将目录拷到服务器后自行执行 pipe.sh
    pipe_sh = phase2_run_dir / "pipe.sh"
    if pipe_sh.is_file():
        try:
            pipe_sh.chmod(pipe_sh.stat().st_mode | 0o111)
        except OSError:
            pass
        print(
            "\n[Phase2] 结构已就绪，未自动运行模拟。请在服务器上进入运行目录后执行：\n"
            f"  cd \"{phase2_run_dir}\"\n"
            "  chmod +x pipe.sh   # 若仍 Permission denied\n"
            "  ./pipe.sh\n"
            "  # 或: bash pipe.sh"
        )

    print("\n[Phase2] Done.")


if __name__ == "__main__":
    main()

