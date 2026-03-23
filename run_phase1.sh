#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# conda 环境 lipid 内：建议用 `python`（而不是 python3），避免 yaml 包缺失
source ~/miniconda3/etc/profile.d/conda.sh
conda activate lipid

python run_phase1_pipeline.py

