# PEG2LNB

为脂质纳气泡添加 PEG 修饰（流程说明见 [`README_zh.md`](README_zh.md)）。

<p align="center">
  <img src="figure.png" alt="PEG2LNB: PEGylated lipid nano-bubble (illustrative rendering)" width="85%">
</p>

Two-phase pipeline for **Martini 3** lipid nano-bubble (LNB) systems with **PEGylated DSPE**, vacuum stripping, re-solvation, bubble-interior solvent removal, and ion neutralization—wired for **GROMACS** (`gmx`).

- **Phase 1**: Build LNB → grow PEG on DSPE → remove water and ions (vacuum-ready system).
- **Phase 2**: Resize box → `gmx solvate` → remove solvent inside the bubble radius (default: `W`) → `gmx genion` to neutralize.

---

## Requirements

| Component | Notes |
|-----------|--------|
| **Python** | 3.10+ recommended |
| **GROMACS** | `gmx` on `PATH` (scripts call `gmx`, not `gmx_mpi`) |
| **Shell** | `csh` or `tcsh` (default `pipe.sh` files are **csh** scripts) |
| **Python packages** | `numpy`, `pyyaml`, `MDAnalysis` |

Optional: **conda** (recommended to isolate dependencies).

---

## Installation

### 1. Clone this repository

```bash
git clone https://github.com/<YOUR_USERNAME>/PEG2LNB.git
cd PEG2LNB
```

### 2. Create a conda environment (recommended)

```bash
conda create -n lipid python=3.11 -y
conda activate lipid
pip install numpy pyyaml MDAnalysis
```

Install **GROMACS** in the same environment or system-wide so that `gmx` is available:

```bash
# Example (conda-forge); adjust for your platform
conda install -c conda-forge gromacs -y
```

On **macOS**, ensure **csh** is available (usually `/bin/csh`).

### 3. Configure paths

Edit the YAML files under `configs/`:

- **`configs/phase1.yaml`** — `lnb_build.*`, `peg.*`, especially `peg.vacuum_run_dir` (where Phase 1 writes `system.gro`, `system.top`, `pipe.sh`, `mdps/`, etc.). The repository defaults to **`output/phase1`** (relative to the repo root).
- **`configs/phase2.yaml`** — `phase1_run_dir`, `phase2_run_dir`, `target_box_nm`, `bubble_radius_nm`, `dry_gro_candidates`, etc. Defaults: **`output/phase1`** and **`output/phase2`**. `dry_gro_candidates` picks the first existing file (e.g. `system_fix.gro` after overlap fixing), otherwise falls back to `system.gro`.

Use **relative paths** under the `PEG2LNB` root when sharing the project publicly (avoid machine-specific absolute paths). The `output/` directory is gitignored so large coordinate files are not committed.

---

## Usage

Run all commands from the `PEG2LNB` directory (or use the provided wrappers).

### Option A — Full pipeline (Phase 1 + Phase 2, optional MD)

`run_all_pipeline.py` runs Phase 1, then Phase 2. By default it also runs `pipe.sh` in each run directory (calls GROMACS). Use `--no-md` to **only** build structures.

```bash
conda activate lipid
python run_all_pipeline.py              # build + run MD in phase1 dir, then phase2 dir
python run_all_pipeline.py --no-md      # build only (no gmx)
python run_all_pipeline.py --skip-phase1   # only Phase 2 (+ optional MD)
python run_all_pipeline.py --skip-phase2   # only Phase 1 (+ optional MD)
```

Convenience wrapper (adjust conda path to match your machine; the script assumes `~/miniconda3`):

```bash
chmod +x run_all.sh
./run_all.sh
./run_all.sh --no-md
```

### Option B — Phase 1 only

```bash
./run_phase1.sh
# or: python run_phase1_pipeline.py
```

Outputs go to `peg.vacuum_run_dir` in `configs/phase1.yaml`.

### Option C — Phase 2 only

Requires a finished Phase 1 directory; paths must match `configs/phase2.yaml`.

```bash
./run_phase2.sh
# or: python run_phase2_pipeline.py
```

### Running MD manually

After each phase, enter the corresponding **run directory** and execute:

```bash
cd /path/to/your/vacuum_or_phase2_run_dir
csh pipe.sh
# or: tcsh pipe.sh
```

(`pipe.sh` is written for **csh**; using `bash pipe.sh` may fail.)

---

## Project layout (high level)

| Path | Role |
|------|------|
| `configs/` | `phase1.yaml`, `phase2.yaml` |
| `default/shared_assets/` | Shared Martini `itps/`, `martini_W.gro` |
| `default/phase1/` | Phase 1 `pipe.sh`, `mdps/` |
| `default/phase2/` | Phase 2 `pipe.sh`, `mdps/` (includes production) |
| `run_*_pipeline.py` | Python drivers for each phase |
| `pipe.py`, `lnb_gener_martini3.py`, … | Core building blocks |

---

## Publishing to GitHub

1. **Scrub secrets and local paths** in `configs/*.yaml` before the first public push.
2. Initialize or use existing git in this folder:

   ```bash
   cd PEG2LNB
   git add README.md README_zh.md .gitignore
   git add default/ configs/ *.py *.sh pipe.py ndx.py ...
   git status   # review carefully
   git commit -m "Initial public release of PEG2LNB"
   ```

3. **Create** an empty repository on GitHub (no README, to avoid conflicts), then:

   ```bash
   git remote add origin https://github.com/<YOUR_USERNAME>/PEG2LNB.git
   git branch -M main
   git push -u origin main
   ```

   If you use the [GitHub CLI](https://cli.github.com/):

   ```bash
   gh auth login
   gh repo create PEG2LNB --public --source=. --remote=origin --push
   ```

---

## License

Add a `LICENSE` file if you plan to redistribute; this repository does not ship one by default.

---

## Citation

If you use this workflow in a publication, cite **GROMACS**, **Martini 3**, and **MDAnalysis** as appropriate, plus your own parameter choices.
