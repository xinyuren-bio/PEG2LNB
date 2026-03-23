# PEG2LNB：LNB-MDT 构建 PEG 纳气泡（去水去离子）-> 再溶剂化加离子（剔除气泡内溶剂）整套流程

该文件夹会把本仓库中与流程相关的脚本与参数“集中”到一起：

- Phase1：LNB-MDT 构建脂质纳气泡 -> DSPE 生长 PEG -> 去除水与溶剂离子（真空体系），**不自动跑 MD**；把运行目录拷到服务器后自行执行 `pipe.sh`
- Phase2：改盒子 -> `gmx solvate` + `gmx genion` -> MDAnalysis 删除气泡质心半径内 `W/NA/CL`，**不自动跑 MD**；同样在服务器上自行执行 `pipe.sh`

## 0. 前置准备

1. 确保你的 conda 环境名为 `lipid`，并且其中包含：
   - `gmx`（或 `gmx_mpi`，脚本按 `gmx` 调用）
   - Python 包：`MDAnalysis`, `numpy`, `pyyaml`

2. 推荐直接用现有命令行：
   - 使用 `python` 而不是 `python3`（你之前遇到过 `python3` 没有 `yaml` 的情况）

## 1. Phase1：构建与去水去离子（不含自动模拟）

默认脚本会读取：
- `PEG2LNB/configs/phase1.yaml`

执行：
```bash
cd /Users/renxinyu/lipid_simulation/PEG2LNB
./run_phase1.sh
```

Phase1 会在 `configs/phase1.yaml` 里 `peg.vacuum_run_dir` 指定的目录（默认示例为 `output/phase1/`，也可改为如 `peg10/`）生成：

- 真空体系 `system.gro` / `system.top`（由根目录 `pipe.py` 产生）
- `mdps/`、`itps/`、`system.ndx`、`pipe.sh`

**在服务器上跑模拟：** 进入该目录后执行（无执行权限时用 `bash pipe.sh`）：

```bash
cd <你的 vacuum_run_dir 绝对路径>
chmod +x pipe.sh   # 可选
./pipe.sh
```

## 2. Phase2：改盒子 + 溶剂化加离子 + 剔除气泡内 W/NA/CL（不含自动模拟）

默认脚本会读取：
- `PEG2LNB/configs/phase2.yaml`

执行：
```bash
cd /Users/renxinyu/lipid_simulation/PEG2LNB
./run_phase2.sh
```

Phase2 会在：
- `PEG2LNB/output/phase2/`

生成：
- `system.gro`（已溶剂化、且剔除了气泡质心半径内的所有 `W/NA/CL`）
- `system.top`、`system.ndx`、`mdps/`、`itps/`、`pipe.sh`

**在服务器上跑模拟：** 进入 `phase2_run_dir`（见 `configs/phase2.yaml`）后同样执行 `./pipe.sh` 或 `bash pipe.sh`。

## 3. 你最可能需要改的参数

1. `PEG2LNB/configs/phase1.yaml` 里：
   - `lnb_build.lipids`（DPPC/DSPE/DOPS/CHOL 数量）
   - `lnb_build.radius / box / area_per_lipid`
   - `peg.peg_length / peg.peg_resname`
   - `lnb_build.gas_num`（气泡内 O2 分子个数）

2. `PEG2LNB/configs/phase2.yaml` 里：
   - `target_box_nm`：修复你遇到的“改了 box 但坐标没平移”的问题（会在溶剂化前对输入体系做平移+wrap）
   - `bubble_radius_nm`：用于删除气泡质心半径内 `W/NA/CL`

## 4. 如果 Phase2 报组选择错误

`gmx genion` 需要你选择要被替换的“水组”。在脚本里默认填的是 `W`（Martini 单珠水残基名）。

如果你的体系水组名不叫 `W`，请在 `phase2.yaml` 中把 `genion_solvent_group` 改成实际提示的组名/组号。

