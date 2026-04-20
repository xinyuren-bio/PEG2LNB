# PEG2LNB：LNB-MDT 构建 PEG 纳气泡（去水去离子）-> 再溶剂化加离子（剔除气泡内溶剂）整套流程

该文件夹会把本仓库中与流程相关的脚本与参数“集中”到一起：

- Phase1：LNB-MDT 构建脂质纳气泡 -> DSPE 生长 PEG -> 去除水与溶剂离子（真空体系），**不自动跑 MD**；把运行目录拷到服务器后自行执行 `pipe.sh`
- Phase2：改盒子 -> `gmx solvate` 加水 -> MDAnalysis 删除气泡质心半径内溶剂（默认仅 `W`）-> `gmx genion -neutral` 置换 `W` 为 `NA`/`CL` 至净电荷为 0（**不按 0.15 M 加盐**），**不自动跑 MD**；同样在服务器上自行执行 `pipe.sh`

## 0. 前置准备

1. 确保你的 conda 环境名为 `lipid`，并且其中包含：
   - `gmx`（或 `gmx_mpi`，脚本按 `gmx` 调用）
   - Python 包：`MDAnalysis`, `numpy`, `pyyaml`

2. 推荐直接用现有命令行：
   - 使用 `python` 而不是 `python3`（你之前遇到过 `python3` 没有 `yaml` 的情况）

3. 默认运行目录（`configs/*.yaml`）已设为相对路径：`output/phase1`、`output/phase2`（首次运行会自动创建；`output/` 已在 `.gitignore` 中忽略，不会提交大块坐标文件）。

## 一键运行（Phase1 + Phase2，可选自动模拟）

```bash
cd PEG2LNB
conda activate lipid
python run_all_pipeline.py          # 默认：两阶段构建后依次在各阶段目录执行 pipe.sh（GROMACS）
python run_all_pipeline.py --no-md  # 仅构建，不跑模拟
```

或使用 `./run_all.sh`（需将脚本内 conda 路径改为你本机的 miniconda/anaconda 路径）。

## 1. Phase1：构建与去水去离子（不含自动模拟）

默认脚本会读取：
- `PEG2LNB/configs/phase1.yaml`

执行：
```bash
cd PEG2LNB
./run_phase1.sh
```

Phase1 会在 `configs/phase1.yaml` 里 `peg.vacuum_run_dir` 指定的目录（默认示例为 `output/phase1/`，也可改为如 `peg10/`）生成：

- 真空体系 `system.gro` / `system.top`（由根目录 `pipe.py` 产生）
- 最终重叠修复后的 `system_fix.gro`（phase1 模拟优先从这个文件开始）
- `mdps/`、`itps/`、`system.ndx`、`pipe.sh`

**在服务器上跑模拟：** 进入该目录后执行（`pipe.sh` 为 **csh** 脚本，请用 `csh`/`tcsh`）：

```bash
cd <你的 vacuum_run_dir>
chmod +x pipe.sh   # 可选
csh pipe.sh
```

## 2. Phase2：改盒子 + 溶剂化 + 剔除气泡内水 + `genion -neutral`（不含自动模拟）

默认脚本会读取：
- `PEG2LNB/configs/phase2.yaml`

执行：
```bash
cd PEG2LNB
./run_phase2.sh
```

Phase2 会在：
- `PEG2LNB/output/phase2/`

生成：
- `system.gro`（已溶剂化、剔除泡内 `W` 后做净电荷中和）
- `system.top`、`system.ndx`、`mdps/`、`itps/`、`pipe.sh`

**在服务器上跑模拟：** 进入 `phase2_run_dir`（见 `configs/phase2.yaml`）后同样执行 `csh pipe.sh`。

## 3. 你最可能需要改的参数

1. `PEG2LNB/configs/phase1.yaml` 里：
   - `lnb_build.lipids`（DPPC/DSPE/DOPS/CHOL 数量）
   - `lnb_build.radius / box / area_per_lipid`
   - `peg.peg_length / peg.peg_resname`
   - `lnb_build.gas_num`（气泡内 O2 分子个数）

2. `PEG2LNB/configs/phase2.yaml` 里：
   - `target_box_nm`：修复你遇到的“改了 box 但坐标没平移”的问题（会在溶剂化前对输入体系做平移+wrap）
   - `bubble_radius_nm` / `remove_resnames_inside_radius`：删泡内溶剂（默认只删 `W`，与「先删水再 genion -neutral」一致）

## 4. 如果 Phase2 报组选择错误

`gmx genion` 需要你选择要被替换的“水组”。在脚本里默认填的是 `W`（Martini 单珠水残基名）。

如果你的体系水组名不叫 `W`，请在 `phase2.yaml` 中把 `genion_solvent_group` 改成实际提示的组名/组号。

## 5. 默认资源目录

当前默认资源已统一放到 `PEG2LNB/default/` 下：

- `default/shared_assets/`：公共 `itps/`、`martini_W.gro`
- `default/phase1/`：phase1 专用 `pipe.sh`、`mdps/`
- `default/phase2/`：phase2 专用 `pipe.sh`、`mdps/`

