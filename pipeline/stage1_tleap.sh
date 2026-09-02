#!/bin/bash
# ============================================================================
# stage1_tleap.sh — 用 tleap 构建 AMBER 拓扑（parm7 + rst7）
#   ff19SB/OL15/TIP3P + 金属 frcmod(12-6) + 12-6-4 离子 + 自动中和 + 12 Å 溶剂化
# 产物: $TOP_DIR/{v.parm7, v.rst7, leap.log}
# ============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
stage_header 1 "tleap 拓扑构建"
init_workdirs
detect_amber
resolve_metal

# 力场名 → leaprc 名
case "$PROT_FF" in
  ff19SB) LP_PROT="leaprc.protein.ff19SB" ;;
  ff14SB) LP_PROT="leaprc.protein.ff14SB" ;;
  *) die "不支持蛋白力场: $PROT_FF（ff19SB / ff14SB）" ;;
esac
case "$DNA_FF" in
  OL15) LP_DNA="leaprc.DNA.OL15" ;;
  bsc1) LP_DNA="leaprc.DNA.bsc1" ;;
  *) die "不支持核酸力场: $DNA_FF（OL15 / bsc1）" ;;
esac
LP_WAT="leaprc.water.$(echo "$WAT_FF" | tr 'A-Z' 'a-z')"
IONS_FRCMOD="${IONS_FRCMOD:-$AMBERHOME/dat/leap/parm/frcmod.ions1lm_1264_tip3p}"
[ -f "$IONS_FRCMOD" ] || die "离子参数不存在: $IONS_FRCMOD（检查 AMBERHOME）"

METAL_FRCMOD_LINE=""
LEAP_PDB="$INPUT_PDB"
if [ "$METAL" = "yes" ]; then
  METAL_FRCMOD_LINE="loadamberparams $FRCMOD"
  # stage0 的 fix_his.py 输出修正副本（HIS→HID/HIE），有则优先用
  [ -f "$TOP_DIR/input_hid.pdb" ] && LEAP_PDB="$TOP_DIR/input_hid.pdb"
fi
# 原子数检查（提前发现 PDB 问题）
NATOM_IN="$(grep -cE '^(ATOM|HETATM)' "$LEAP_PDB")"
info "输入 PDB 原子数: $NATOM_IN"

banner "渲染 leap.in → $TOP_DIR/leap.in"
render "$TPL_DIR/leap.in.tmpl" "$TOP_DIR/leap.in" \
  LP_PROT="$LP_PROT" LP_DNA="$LP_DNA" LP_WAT="$LP_WAT" \
  IONS_FRCMOD="$IONS_FRCMOD" METAL_FRCMOD_LINE="$METAL_FRCMOD_LINE" \
  INPUT_PDB="$(cd "$RUN_DIR" && realpath "$LEAP_PDB" 2>/dev/null || echo "$LEAP_PDB")" \
  BOX_DIST="$BOX_DIST" OUT_PARM7="$TOP_DIR/v.parm7" OUT_RST7="$TOP_DIR/v.rst7"

info "运行 tleap（约 1-3 分钟）"
exec_cmd "'$TLEAP' -f '$TOP_DIR/leap.in' > '$TOP_DIR/leap.log' 2>&1"
grep -q "Error" "$TOP_DIR/leap.log" && { err "tleap 报错（见 $TOP_DIR/leap.log 末尾）"; tail -20 "$TOP_DIR/leap.log" >&2; exit 1; }
[ -f "$TOP_DIR/v.parm7" ] && [ -f "$TOP_DIR/v.rst7" ] || die "tleap 未产出 parm7/rst7（看日志）"

# 校验: Mn 电荷 = +2 / 总电荷 ≈ 0 / 原子数
"$PYTHON" - "$TOP_DIR/v.parm7" <<'PYEOF' || die "拓扑电荷校验失败"
import sys
import parmed as pmd
t = pmd.load_file(sys.argv[1])
tot = sum(a.charge for a in t.atoms)
mn = [a for a in t.atoms if a.name == "MN"]
qmn = sum(a.charge for a in mn) if mn else None
print(f"  n_atom={len(t.atoms)}  total_charge={tot:+.2f}"
      + (f"  Mn_charge={qmn:+.2f}" if qmn is not None else ""))
if mn and abs(qmn - 2.0) > 0.01:
    sys.exit(1)
if abs(tot) > 0.5:
    print("  [WARN] 总电荷偏离 0 —— 检查离子中和", file=sys.stderr)
PYEOF
ok "拓扑构建完成: $TOP_DIR/v.parm7 (+v.rst7)"
