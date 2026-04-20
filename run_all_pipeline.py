# ==================================================
# 功能说明：串联 Phase1 与 Phase2 构建流程，可选在各阶段结束后自动执行 pipe.sh 运行 GROMACS。
# 使用方法：在 PEG2LNB 目录下执行 `python run_all_pipeline.py`；加 `--no-md` 则只构建不跑模拟；`--skip-phase1` / `--skip-phase2` 可跳过对应阶段。
# 依赖环境：与 run_phase1_pipeline.py / run_phase2_pipeline.py 相同（Python、PyYAML、MDAnalysis、numpy、gmx）；执行模拟需系统中有 csh 或 tcsh。
# 生成时间：2026-04-20
# ==================================================

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def _load_yaml_min(path: Path) -> dict:
    try:
        import yaml  # type: ignore

        return yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except Exception:
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


def _resolve_phase1_run_dir(repo_root: Path, cfg: dict) -> Path:
    peg = cfg.get("peg") or {}
    rel = peg.get("vacuum_run_dir", "")
    p = Path(rel)
    return p if p.is_absolute() else (repo_root / p).resolve()


def _resolve_phase2_run_dir(repo_root: Path, cfg: dict) -> Path:
    rel = cfg.get("phase2_run_dir", "")
    p = Path(rel)
    return p if p.is_absolute() else (repo_root / p).resolve()


def _run_pipe_sh(run_dir: Path) -> None:
    pipe_sh = run_dir / "pipe.sh"
    if not pipe_sh.is_file():
        raise FileNotFoundError(f"未找到 {pipe_sh}，无法运行模拟。")
    csh = shutil.which("csh") or shutil.which("tcsh")
    if not csh:
        raise RuntimeError("未在 PATH 中找到 csh 或 tcsh；pipe.sh 为 csh 脚本，请安装后再试。")
    print(f"\n[run_all] 在目录中执行模拟: {run_dir}\n  {csh} pipe.sh\n")
    subprocess.run([csh, "pipe.sh"], cwd=str(run_dir), check=True)


def main() -> None:
    repo_root = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(
        description="PEG2LNB：Phase1 + Phase2 一键流水线；默认在每阶段构建完成后执行该阶段 pipe.sh。"
    )
    parser.add_argument(
        "--no-md",
        action="store_true",
        help="仅运行 Python 构建步骤，不执行各阶段目录下的 pipe.sh（GROMACS）。",
    )
    parser.add_argument(
        "--skip-phase1",
        action="store_true",
        help="跳过 Phase1（假定 vacuum_run_dir 中已有 Phase1 产物）。",
    )
    parser.add_argument(
        "--skip-phase2",
        action="store_true",
        help="跳过 Phase2。",
    )
    args = parser.parse_args()
    run_md = not args.no_md

    cfg1_path = repo_root / "configs" / "phase1.yaml"
    cfg2_path = repo_root / "configs" / "phase2.yaml"
    cfg1 = _load_yaml_min(cfg1_path)
    cfg2 = _load_yaml_min(cfg2_path)
    phase1_run_dir = _resolve_phase1_run_dir(repo_root, cfg1)
    phase2_run_dir = _resolve_phase2_run_dir(repo_root, cfg2)

    if not args.skip_phase1:
        print("\n[run_all] >>> Phase1 构建 …\n")
        subprocess.run(
            [sys.executable, str(repo_root / "run_phase1_pipeline.py")],
            check=True,
            cwd=str(repo_root),
        )
        if run_md:
            _run_pipe_sh(phase1_run_dir)
    else:
        print("\n[run_all] 已跳过 Phase1 构建。")

    if not args.skip_phase2:
        print("\n[run_all] >>> Phase2 构建 …\n")
        subprocess.run(
            [sys.executable, str(repo_root / "run_phase2_pipeline.py")],
            check=True,
            cwd=str(repo_root),
        )
        if run_md:
            _run_pipe_sh(phase2_run_dir)
    else:
        print("\n[run_all] 已跳过 Phase2 构建。")

    print("\n[run_all] 全部完成。")


if __name__ == "__main__":
    main()
