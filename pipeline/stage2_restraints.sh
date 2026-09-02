#!/bin/bash
# ============================================================================
# stage2_restraints.sh — 生成约束 mask + DISANG + 金属配位体检
#   make_restraints.py 读 parm7/rst7:
#     - 自动识别 蛋白/DNA/Mn 残基区间 → masks.env (HEAVY_BB_MASK/BB_MASK/PROT_RANGE/DNA_RANGE/MN_RES)
#     - 自动探测 Mn 3.5 Å 内 His/Glu 配位原子 → disang.txt flat-bottom 约束
#     - 输出 COORD1..4（配位原子 res:atom，stage3/4 用）
# 金属体系额外用 cpptraj 验证 rst7 中配位距离 ≈ 1.9-2.6 Å（三齿健康）
# ============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
stage_header 2 "约束与金属配位体检"
init_workdirs
detect_amber
resolve_metal
[ -f "$TOP_DIR/v.parm7" ] || die "先跑阶段 1（缺 $TOP_DIR/v.parm7）"

info "运行 make_restraints.py（自动探测掩码与配位）"
exec_cmd "'$PYTHON' '$TOOLS_DIR/make_restraints.py' '$TOP_DIR/v.parm7' '$TOP_DIR/v.rst7' '$TOP_DIR'"
[ -f "$TOP_DIR/masks.env" ] || die "masks.env 未生成"
source "$TOP_DIR/masks.env"

# 无金属体系：写一个空 DISANG 占位（Amber 要求 DISANG 路径存在但 nmropt=0 不读）
if [ "$METAL" = "no" ]; then
  : > "$TOP_DIR/disang.txt"
  info "无金属：DISANG 置空（阶段 3 模板将不启用 nmropt）"
fi

if [ "$METAL" = "yes" ]; then
  [ "${MN_RES:-}" = "-" ] && die "METAL=yes 但 masks.env MN_RES=-（未识别到 Mn 残基）"
  info "Mn 残基=$MN_RES  配位原子: ${COORD1:-?} ${COORD2:-} ${COORD3:-} ${COORD4:-}"
  # cpptraj 体检：rst7 中 4 条配位距离
  cat > "$TOP_DIR/mn_check.in" <<EOF
parm $TOP_DIR/v.parm7
trajin $TOP_DIR/v.rst7
distance c1 :${MN_RES}@MN :${COORD1%%:*}@${COORD1##*:} out $TOP_DIR/mn_c1.dat
distance c2 :${MN_RES}@MN :${COORD2%%:*}@${COORD2##*:} out $TOP_DIR/mn_c2.dat
distance c3 :${MN_RES}@MN :${COORD3%%:*}@${COORD3##*:} out $TOP_DIR/mn_c3.dat
distance c4 :${MN_RES}@MN :${COORD4%%:*}@${COORD4##*:} out $TOP_DIR/mn_c4.dat
run
quit
EOF
  exec_cmd "'$CPPTRAJ' -i '$TOP_DIR/mn_check.in' > '$TOP_DIR/mn_check.log' 2>&1"
  banner "Mn 配位距离体检（rst7 起始结构）"
  for c in c1 c2 c3 c4; do
    if [ -f "$TOP_DIR/mn_$c.dat" ]; then
      d="$(tail -1 "$TOP_DIR/mn_$c.dat" | awk '{print $2}')"
      ok "$c = ${d} Å"
    fi
  done
  warn "以上距离应落在 1.9-2.6 Å（三齿健康）。若 >3.5 Å：输入 PDB 中 Mn 未与 HEXXH 配位，"
  warn "请检查/重建起始结构再重跑（详见 docs/METAL_CHECKLIST.md）。"
fi
ok "约束生成完毕: masks.env + disang.txt"
