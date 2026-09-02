#!/bin/bash
# ============================================================================
# stage5_mmpbsa.sh — MM-PBSA/GBSA 结合自由能 + per-residue 分解（可选，默认关）
#   DO_MMPBSA=yes 才执行；产物拷入 results/$SYSTEM/ 供 fig07/fig21/fig24 使用
#   体系要求: 蛋白 + DNA 复合物（rec=蛋白, lig=DNA, cp=复合物）。Mn 不参与
#     GB 计算（sander 不识别 Mn 原子类型，原项目铁律）→ rec/cp 均剥 Mn。
#   方法: MMPBSA.py GB(igb=5, saltcon=0.15) 做 per-residue 分解 + PB(istrng=0.15) 总自由能
#   帧采样: 每 rep 默认取 1..last, interval=10（1500 帧 → 150 帧，文献常见下限）
#   并发: MMPBSA.py 对临时文件敏感 → 每 rep 单独 cd 后执行（历史教训）
# ============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
stage_header 5 "MM-PBSA (可选)"
init_workdirs
detect_amber
resolve_metal

if [ "${DO_MMPBSA:-no}" != "yes" ]; then
  info "DO_MMPBSA=no — 跳过（fig07/fig21/fig24 将缺 MM-PBSA 数据）"
  exit 0
fi
[ -f "$TOP_DIR/v.parm7" ] || die "先跑阶段 1"
[ -f "$TOP_DIR/masks.env" ] || die "先跑阶段 2"
source "$TOP_DIR/masks.env"

# 体系门槛：需同时有蛋白与 DNA
if [ "${PROT_RANGE:-}" = "-" ] || [ "${DNA_RANGE:-}" = "-" ]; then
  warn "MM-PBSA 需要 蛋白+DNA 受体-配体对（当前 PROT=$PROT_RANGE DNA=$DNA_RANGE）— 跳过"
  exit 0
fi

PROT_MASK=":$PROT_RANGE"     # e.g. :1-254
DNA_MASK=":$DNA_RANGE"       # e.g. :255-278
MN_MASK=":MN"
[ "${MN_RES:-}" != "-" ] && MN_MASK=":$MN_RES" || MN_MASK=""

MMPBSA_BIN="${MMPBSA_BIN:-$AMBERHOME/bin/MMPBSA.py}"
[ -x "$MMPBSA_BIN" ] || die "未找到 MMPBSA.py（$MMPBSA_BIN）—— 需 AmberTools 完整安装"
MMPBSA_INTERVAL="${MMPBSA_INTERVAL:-10}"

# 用 parmed 从 v.parm7 剥出 rec/lig/cp（坐标系一致: strip 后仍写 rst7 供参考）
banner "准备 rec(蛋白)/lig(DNA)/cp(复合物) 拓扑（parmed strip）"
for kind in rec lig cp; do
  [ -f "$TOP_DIR/$kind.parm7" ] || {
    cat > "$TOP_DIR/$kind.in" <<EOF
parm $TOP_DIR/v.parm7
reference $TOP_DIR/v.rst7
EOF
    # strip 序列：顺序 = 反向剥离（后写的先 strip 无妨，parmed 允许任意序）
    case "$kind" in
      rec) echo "strip $DNA_MASK" >> "$TOP_DIR/$kind.in" ;;
      lig) echo "strip $PROT_MASK" >> "$TOP_DIR/$kind.in" ;;
      cp)  : ;;
    esac
    # 溶剂离子与 Mn 一律剥除（cp 也剥 Mn → GB 不识别 Mn）
    echo "strip :WAT,:Na+,:Cl-" >> "$TOP_DIR/$kind.in"
    [ -n "$MN_MASK" ] && echo "strip $MN_MASK" >> "$TOP_DIR/$kind.in"
    cat >> "$TOP_DIR/$kind.in" <<EOF
outparm $TOP_DIR/$kind.parm7 $TOP_DIR/$kind.rst7
quit
EOF
    info "parmed → $kind.parm7"
    exec_cmd "'$AMBERHOME/bin/parmed' -i '$TOP_DIR/$kind.in' > '$TOP_DIR/$kind.parmed.log' 2>&1" \
      || { tail -10 "$TOP_DIR/$kind.parmed.log" >&2; die "$kind strip 失败"; }
  }
done

# 逐 rep 跑 MMPBSA.py
for rep in $(seq 1 "$NREP"); do
  W="$OUT_DIR/rep$rep"
  [ -s "$W/prod.nc" ] || { info "rep$rep 无 prod.nc — 跳过"; continue; }
  [ -s "$W/FINAL_RESULTS.dat" ] && { ok "rep$rep MMPBSA 已有产物 — 跳过"; continue; }

  # 用 python 数 prod.nc 帧数（cpptraj 亦可；轻量 python 直读 netcdf 头不划算，
  # 直接信任 nstlim/ntwx: 帧数 = PROD_NS*10 对应 100ps/帧，见 prod 模板注释）
  END_FRAME=$((PROD_NS * 10))
  render "$TPL_DIR/mmpbsa.in.tmpl" "$W/mmpbsa.in" \
    START=1 END="$END_FRAME" INTERVAL="$MMPBSA_INTERVAL"
  banner "rep$rep MM-PBSA (frames 1-${END_FRAME} step ${MMPBSA_INTERVAL})"
  info "运行 MMPBSA.py（GB per-residue + PB 总量）—— 可能数小时"
  # 历史教训: MMPBSA.py 临时文件随 CWD 生成 → cd $W 防并发污染
  ( cd "$W" && exec_cmd "'$MMPBSA_BIN' -O -i '$W/mmpbsa.in' \
      -o '$W/FINAL_RESULTS.dat' -do '$W/FINAL_DECOMP_MMPBSA.dat' \
      -sp '$TOP_DIR/v.parm7' \
      -cp '$TOP_DIR/cp.parm7' -rp '$TOP_DIR/rec.parm7' -lp '$TOP_DIR/lig.parm7' \
      -y '$W/prod.nc' > '$W/mmpbsa.log' 2>&1" ) \
    || { tail -20 "$W/mmpbsa.log" >&2; die "rep$rep MMPBSA 失败"; }
  ok "rep$rep ΔG 完成"
done

# rep1 结果拷入 results/ 顶层（fig 脚本硬编码读取位置）
R1="$OUT_DIR/rep1"
if [ -s "$R1/FINAL_RESULTS.dat" ]; then
  cp -f "$R1/FINAL_RESULTS.dat" "$RES_DIR/FINAL_RESULTS.dat"
  cp -f "$R1/FINAL_DECOMP_MMPBSA.dat" "$RES_DIR/FINAL_DECOMP_MMPBSA.dat" 2>/dev/null || \
    warn "rep1 无 per-residue 分解（fig21/fig24 将缺）"
  ok "结果已汇总 → $RES_DIR/FINAL_RESULTS.dat"
else
  warn "rep1 无 MMPBSA 产物 — results/ 顶层未汇总"
fi
