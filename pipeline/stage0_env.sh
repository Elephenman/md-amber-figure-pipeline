#!/bin/bash
# ============================================================================
# stage0_env.sh — 环境体检 + 输入 PDB 体检
#   1) Amber/pmemd/cpptraj/python 探测     2) 输入 PDB 存在/链/金属检查
#   3) 金属体系: HIS→HID 质子化态改写（Mn 配位必需）
# 失败即退出；本阶段不产生 .done（每次重跑都重新体检）
# ============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
stage_header 0 "环境与输入体检"
init_workdirs

info "1/3 引擎探测"
detect_amber
need "$PYTHON"
"$PYTHON" -c "import numpy, scipy, matplotlib" 2>/dev/null \
  || warn "python3 缺少 numpy/scipy/matplotlib（出图阶段需要；跑 MD 可忽略）"
nvidia-smi >/dev/null 2>&1 && ok "检测到 NVIDIA GPU" || info "未检测到 GPU（将用 CPU 引擎）"

info "2/3 输入 PDB 体检"
[ -f "$INPUT_PDB" ] || die "输入 PDB 不存在: $INPUT_PDB"
resolve_metal
"$PYTHON" "$TOOLS_DIR/check_input_pdb.py" "$INPUT_PDB" --metal "$METAL"

info "3/3 金属配位 HIS → HID/HIE（ε/δ 去质子化，按 Mn 配位原子自动判定）"
if [ "$METAL" = "yes" ]; then
  # 输出到 TOP_DIR 修正副本（不污染用户输入）；stage1 自动优先用 input_hid.pdb
  "$PYTHON" "$TOOLS_DIR/fix_his.py" "$INPUT_PDB" "$TOP_DIR/input_hid.pdb"
else
  info "无金属体系：跳过 HID 改写"
fi
ok "体检通过"
