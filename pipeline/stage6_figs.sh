#!/bin/bash
# ============================================================================
# stage6_figs.sh — 一键 30 图（fig01-30）+ HTML 画廊
#   前置: stage4 数据 (results/$SYSTEM/*.dat + cpptraj_raw/)、MM-PBSA(可选 fig07/21/24)
#   调用: python scripts/plot_pretty_figs{1,2,3}.py  +  figures_gallery.py
#   体系门槛: 30 图脚本面向「蛋白+核酸(+金属)」复合物体系（WT·S1·Mn 案例模板）。
#     无 DNA（apo/freedna/裸蛋白）→ 整体跳过出图（fig01-08 也含 DNA 视图）。
#     plot2/3 (fig09-30) 为 PprI 案例深度绑定脚本（残基/序列注释）；
#     复现 PprI 系体系可直接用，新蛋白体系建议 DO_FIGS=no 或仅取 plot1 通用图。
#   数据预检: 缺关键文件必报错（缺哪个、补跑哪个 stage 都会提示）；
#     可选数据（bfactor_prot.dat / MM-PBSA decomp）已在脚本内做占位守卫。
# ============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
stage_header 6 "一键 30 图 + 画廊"
init_workdirs

if [ "${DO_FIGS:-yes}" != "yes" ]; then
  info "DO_FIGS=no — 跳过出图"
  exit 0
fi

# python + 绘图库检查
need "$PYTHON"
"$PYTHON" -c "import numpy, scipy, matplotlib" 2>/dev/null \
  || die "缺少 numpy/scipy/matplotlib（先装: pip install numpy scipy matplotlib）"

RAW="$RES_DIR/cpptraj_raw"
[ -d "$RAW" ] || die "缺 $RAW —— 先跑阶段 4"
PROT_RANGE_SRC="$TOP_DIR/masks.env"
[ -f "$PROT_RANGE_SRC" ] || die "缺 masks.env —— 先跑阶段 2"
source "$PROT_RANGE_SRC"
export MDEASY_RES="$RES_DIR"
export MDEASY_CASE="$SYSTEM"

# ---------- 体系类型门控 ----------
# 图集(30 张)面向「蛋白+核酸(+金属)」复合物; apo/freedna(无 DNA) 体系
# 连 fig01-08 也含 DNA 视图(rmsd_dna/rg_dna/sasa_dna), 故整体不适配 → 跳过出图。
if [ "${DNA_RANGE:-}" = "-" ] || [ -z "${DNA_RANGE:-}" ]; then
  info "本体系无 DNA（apo/freedna/裸蛋白）——30 图脚本面向蛋白+核酸复合物，跳过出图"
  info "  （裸蛋白体系如需二级结构图, 可手工跑 cpptraj dssp; 或用 DO_FIGS=no 跳过）"
  exit 0
fi

# ---------- 必需数据预检（stage4 顶层 .dat + cpptraj_raw 关键量） ----------
banner "数据预检"
REQUIRED_TOP="rmsd_prot.dat rmsf_prot.dat rg_prot.dat sasa_prot.dat sasa_prot_iso.dat sasa_cplx.dat pca_proj.dat pca_all.dat prod.out"
REQUIRED_RAW="phi.dat psi.dat dssp.dat dssp_sum.dat"
MISSING=""
for f in $REQUIRED_TOP; do [ -s "$RES_DIR/$f" ] || MISSING="$MISSING $f(top)"; done
for f in $REQUIRED_RAW; do [ -s "$RAW/$f" ] || MISSING="$MISSING $f(raw)"; done

# DNA 分支文件（an2b contacts/hbond + an2d/an4 DNA 形态 + an2c pock）
for f in rmsd_dna.dat rg_dna.dat rmsf_dna.dat sasa_dna.dat sasa_dna_iso3.dat; do
  [ -s "$RES_DIR/$f" ] || MISSING="$MISSING $f(dna)"
done
for f in nat_series.dat nat_res.dat hb_all_avg.dat hb_all_series.dat; do
  [ -s "$RAW/$f" ] || MISSING="$MISSING $f(raw,dna)"
done
[ -s "$RAW/dna_ee.dat" ]  || MISSING="$MISSING dna_ee.dat(raw,dna)"
[ -s "$RAW/dna_bend.dat" ] || MISSING="$MISSING dna_bend.dat(raw,dna)"
[ -s "$RAW/pock_rmsd.dat" ] || MISSING="$MISSING pock_rmsd.dat(raw)"
# 金属分支（fig08/30 + dccm；bfactor_prot.dat 可选 → fig28 已做占位守卫）
if [ "${MN_RES:-}" != "-" ] && [ -n "${MN_RES:-}" ]; then
  for f in mn_H71.dat mn_H75.dat mn_E102a.dat mn_E102b.dat mn_rdf.dat; do
    [ -s "$RES_DIR/$f" ] || MISSING="$MISSING $f(metal)"
  done
  [ -s "$RAW/dccm_prot.dat" ] || MISSING="$MISSING dccm_prot.dat(raw)"
fi
[ -z "$MISSING" ] || die "缺必需数据:$MISSING\n  → 补跑: bash md_easy.sh -c <config> --from 4   （MM-PBSA 图需阶段 5）"

# ---------- 执行三连绘图 ----------
banner "fig01-08 (plot1)"
( cd "$SCRIPTS_DIR" && "$PYTHON" plot_pretty_figs.py ) || die "plot1 失败（见上方 traceback）"

banner "fig09-22 (plot2)"
( cd "$SCRIPTS_DIR" && "$PYTHON" plot_pretty_figs2.py ) || die "plot2 失败"
banner "fig23-30 (plot3)"
( cd "$SCRIPTS_DIR" && "$PYTHON" plot_pretty_figs3.py ) || die "plot3 失败"

# ---------- 画廊 ----------
banner "HTML 画廊"
if [ -f "$SCRIPTS_DIR/figures_gallery.py" ]; then
  ( cd "$SCRIPTS_DIR" && "$PYTHON" figures_gallery.py ) || warn "画廊生成失败（不影响图）"
fi
FIGDIR="$RES_DIR/figures_pretty"
nfig="$(ls "$FIGDIR"/fig*.png 2>/dev/null | wc -l)"
ok "完成: $nfig 张图 → $FIGDIR/  （打开 index.html 看画廊）"
