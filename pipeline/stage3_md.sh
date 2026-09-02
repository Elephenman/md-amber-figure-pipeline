#!/bin/bash
# ============================================================================
# stage3_md.sh — 生产 MD 主流程：min1 → min2 → heat → equil → prod
#   每个重复一个目录 out/$SYSTEM/rep$REP/；GPU 失败自动退 CPU；单阶段断点续跑
#   min1 : 骨架+金属约束最小化(25 kcal/mol·Å², 3000 cyc)     —— 缓解输入 clash
#   min2 : 全体系无约束最小化(5000 cyc)
#   heat : NVT 0→300 K 50000 步(100 ps)，骨架 wt=10 + Mn DISANG
#   equil: NPT 300 K 250000 步(500 ps)，骨架 wt=2 + Mn DISANG
#   prod : NPT 300 K 自由动力学(ntr=0) nstlim=500000×PROD_NS，仅 Mn DISANG
# 模板占位符: 由 common.sh render() 渲染（见 templates/*.in.tmpl）
# ============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
stage_header 3 "MD (min1→min2→heat→equil→prod)"
init_workdirs
detect_amber
resolve_metal
[ -f "$TOP_DIR/v.parm7" ] || die "先跑阶段 1"
[ -f "$TOP_DIR/masks.env" ] || die "先跑阶段 2"
source "$TOP_DIR/masks.env"

# DISANG / nmropt 渲染块（金属体系启用 nmropt + DISANG；无金属整块清空）
NMRPT_LINE=""          # cntrl 内 nmropt 行
DISANG_TAIL=""         # 文件尾 &wt/DISANG 块
if [ "$METAL" = "yes" ]; then
  NMRPT_LINE="  nmropt=1,"
  DISANG_TAIL=$'&wt type=\x27END\x27\n/\nDISANG='"$TOP_DIR/disang.txt"
fi

NSTLIM=$((PROD_NS * 500000))
BB_SELECT="${BB_MASK:-$HEAVY_BB_MASK}"
HB_SELECT="${HEAVY_BB_MASK:-$BB_SELECT}"

for rep in $(seq 1 "$NREP"); do
  W="$OUT_DIR/rep$rep"
  mkdir -p "$W"
  banner "== $SYSTEM rep$rep =="

  # --- 断点：prod.nc 已在 → 跳过该重复 ---
  if [ -s "$W/prod.nc" ] && [ "$(wc -c < "$W/prod.nc")" -gt 1000000 ]; then
    ok "rep$rep 已产出 prod.nc，跳过（删除可重跑）"
    continue
  fi

  # 用 rst7（若该重复已有 equil.rst 则续跑热化后阶段）
  RST="$TOP_DIR/v.rst7"
  for st in min1 min2 heat equil; do
    [ -s "$W/$st.rst" ] && RST="$W/$st.rst"
  done

  # ---------- min1（骨架约束最小化，缓解输入 clash） ----------
  if [ ! -s "$W/min1.rst" ]; then
    render "$TPL_DIR/min1.in.tmpl" "$W/min1.in" RESTRAINT_WT="$MIN_WT" RESTRAINT_MASK="$HB_SELECT"
    info "min1: 骨架约束最小化 (wt=$MIN_WT)"
    exec_cmd "'$MD_ENGINE' -O -i '$W/min1.in' -p '$TOP_DIR/v.parm7' -c '$RST' -o '$W/min1.out' -r '$W/min1.rst' -ref '$RST' || { warn 'GPU min1 失败，退回 CPU'; '$PMEMD_CPU' -O -i '$W/min1.in' -p '$TOP_DIR/v.parm7' -c '$RST' -o '$W/min1.out' -r '$W/min1.rst' -ref '$RST'; }"
    [ -s "$W/min1.rst" ] || die "min1 失败（见 $W/min1.out）"
  else ok "min1 已有产物，跳过"; fi

  # ---------- min2（无约束最小化） ----------
  if [ ! -s "$W/min2.rst" ]; then
    render "$TPL_DIR/min2.in.tmpl" "$W/min2.in"
    info "min2: 无约束最小化"
    exec_cmd "'$MD_ENGINE' -O -i '$W/min2.in' -p '$TOP_DIR/v.parm7' -c '$W/min1.rst' -o '$W/min2.out' -r '$W/min2.rst' || { warn 'GPU min2 失败，退回 CPU'; '$PMEMD_CPU' -O -i '$W/min2.in' -p '$TOP_DIR/v.parm7' -c '$W/min1.rst' -o '$W/min2.out' -r '$W/min2.rst'; }"
    [ -s "$W/min2.rst" ] || die "min2 失败（见 $W/min2.out）"
  else ok "min2 已有产物，跳过"; fi

  # ---------- heat（NVT 0→TEMP，骨架 wt=HEAT_WT + DISANG） ----------
  if [ ! -s "$W/heat.rst" ]; then
    render "$TPL_DIR/heat.in.tmpl" "$W/heat.in" \
      RESTRAINT_WT="$HEAT_WT" RESTRAINT_MASK="$HB_SELECT" \
      NMRPT_LINE="$NMRPT_LINE" DISANG_TAIL="$DISANG_TAIL" TEMP="$TEMP"
    info "heat: NVT 0→${TEMP} K (100 ps, 骨架 wt=$HEAT_WT)"
    exec_cmd "'$MD_ENGINE' -O -i '$W/heat.in' -p '$TOP_DIR/v.parm7' -c '$W/min2.rst' -ref '$W/min2.rst' -o '$W/heat.out' -r '$W/heat.rst' -x '$W/heat.nc' -inf '$W/heat.mdinfo' || { warn 'GPU heat 失败，退回 CPU'; '$PMEMD_CPU' -O -i '$W/heat.in' -p '$TOP_DIR/v.parm7' -c '$W/min2.rst' -ref '$W/min2.rst' -o '$W/heat.out' -r '$W/heat.rst' -x '$W/heat.nc' -inf '$W/heat.mdinfo'; }"
    [ -s "$W/heat.rst" ] || die "heat 失败（见 $W/heat.out）"
  else ok "heat 已有产物，跳过"; fi

  # ---------- equil（NPT 300 K，骨架 wt=EQUIL_WT + DISANG） ----------
  if [ ! -s "$W/equil.rst" ]; then
    render "$TPL_DIR/equil.in.tmpl" "$W/equil.in" \
      RESTRAINT_WT="$EQUIL_WT" RESTRAINT_MASK="$BB_SELECT" \
      NMRPT_LINE="$NMRPT_LINE" DISANG_TAIL="$DISANG_TAIL" TEMP="$TEMP"
    info "equil: NPT ${TEMP} K (500 ps, 骨架 wt=$EQUIL_WT)"
    exec_cmd "'$MD_ENGINE' -O -i '$W/equil.in' -p '$TOP_DIR/v.parm7' -c '$W/heat.rst' -ref '$W/heat.rst' -o '$W/equil.out' -r '$W/equil.rst' -x '$W/equil.nc' -inf '$W/equil.mdinfo' || { warn 'GPU equil 失败，退回 CPU'; '$PMEMD_CPU' -O -i '$W/equil.in' -p '$TOP_DIR/v.parm7' -c '$W/heat.rst' -ref '$W/heat.rst' -o '$W/equil.out' -r '$W/equil.rst' -x '$W/equil.nc' -inf '$W/equil.mdinfo'; }"
    [ -s "$W/equil.rst" ] || die "equil 失败（见 $W/equil.out）"
  else ok "equil 已有产物，跳过"; fi

  # ---------- prod（NPT 自由动力学 + DISANG，nstlim=NSTLIM） ----------
  if [ ! -s "$W/prod.nc" ]; then
    render "$TPL_DIR/prod.in.tmpl" "$W/prod.in" \
      NSTLIM="$NSTLIM" NMRPT_LINE="$NMRPT_LINE" DISANG_TAIL="$DISANG_TAIL" TEMP="$TEMP"
    info "prod: NPT ${TEMP} K 自由动力学 ${PROD_NS} ns (ntr=0, nstlim=$NSTLIM)"
    exec_cmd "'$MD_ENGINE' -O -i '$W/prod.in' -p '$TOP_DIR/v.parm7' -c '$W/equil.rst' -ref '$W/equil.rst' -o '$W/prod.out' -r '$W/prod.rst' -x '$W/prod.nc' -inf '$W/prod.mdinfo' || { warn 'GPU prod 失败，退回 CPU（检查 prod.out 尾部）'; tail -20 '$W/prod.out' >&2 || true; exit 1; }"
    [ -s "$W/prod.nc" ] || die "prod 失败（见 $W/prod.out）"
  else ok "prod 已有产物，跳过"; fi

  # ---------- 快速健康检查：RESTRAINT 能量 & 温度 ----------
  if [ "$METAL" = "yes" ]; then
    RE="$(grep -m1 'RESTRAINT' "$W/prod.out" | awk '{print $3}')"
    if [ -n "$RE" ] && awk -v x="$RE" 'BEGIN{exit !(x>100)}'; then
      err "prod 首帧 RESTRAINT 能量 $RE kcal/mol —— Mn 配位约束崩坏（>100）"
      err "排查: 看阶段 2 配位体检与 docs/METAL_CHECKLIST.md（常见：起始 DNA 贴近 Mn）"
      exit 1
    fi
    [ -n "$RE" ] && ok "prod 首帧 RESTRAINT 能量 = $RE kcal/mol（健康应 <6）"
  fi
  ok "rep$rep 完成 → $W/prod.nc"
done
ok "MD 阶段全部完成"
