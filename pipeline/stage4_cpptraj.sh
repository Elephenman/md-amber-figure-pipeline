#!/bin/bash
# ============================================================================
# stage4_cpptraj.sh — 轨迹数据采集（cpptraj，CPU）
#   输出: results/$SYSTEM/            ← an0 基础量 (rmsd/rg/rmsf/sasa/pca/mn 配位时序)
#         results/$SYSTEM/cpptraj_raw ← an1..an5 专项量 (dssp/contacts/hbond/pock/ee/bend/dccm/rdf)
#         results/$SYSTEM/prod.out    ← rep1 生产日志 (fig09)
# 体系类型自适应: apo(无DNA)/freedna(无蛋白)/无金属 自动跳过无关分析
# ============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
stage_header 4 "cpptraj 轨迹数据采集"
init_workdirs
detect_amber
resolve_metal
[ -f "$TOP_DIR/v.parm7" ] || die "先跑阶段 1"
[ -f "$TOP_DIR/masks.env" ] || die "先跑阶段 2"
source "$TOP_DIR/masks.env"
RAW="$RES_DIR/cpptraj_raw"
mkdir -p "$RAW"

# 找到至少一个 prod.nc 作为主轨迹（多重复全部纳入）
TRAJS=""
for rep in $(seq 1 "$NREP"); do
  [ -s "$OUT_DIR/rep$rep/prod.nc" ] && TRAJS="$TRAJS $OUT_DIR/rep$rep/prod.nc"
done
[ -n "$TRAJS" ] || die "未找到 prod.nc（先跑阶段 3）"
TRJ1="${TRAJS%% *}"
# rep1 的 prod.out 归档给 fig09
cp -f "$OUT_DIR/rep1/prod.out" "$RES_DIR/prod.out" 2>/dev/null || warn "rep1/prod.out 不存在（fig09 将缺数据）"

# 派生几何
PROT_P=""; DNA_P=""; PROT_HAS=0; DNA_HAS=0
if [ "${PROT_RANGE:-}" != "-" ] && [ -n "${PROT_RANGE:-}" ]; then PROT_P=":$PROT_RANGE"; PROT_HAS=1; fi
if [ "${DNA_RANGE:-}" != "-" ] && [ -n "${DNA_RANGE:-}" ]; then DNA_P=":$DNA_RANGE"; DNA_HAS=1; fi
if [ "$PROT_HAS" = "0" ] && [ "$DNA_HAS" = "0" ]; then die "既无蛋白也无 DNA 掩码，无法分析"; fi

TRAJIN_BLOCK=""
for t in $TRAJS; do TRAJIN_BLOCK="$TRAJIN_BLOCK
trajin $t 1 last 1"; done

# Mn 配位原子掩码 (pock RMSD 用)
MN_POCK=""
COORD_RES=""
if [ "$METAL" = "yes" ] && [ "${MN_RES:-}" != "-" ]; then
  for i in 1 2 3 4; do
    c=$(eval "echo \${COORD$i:-}")
    [ "$c" = "-" ] && continue
    r="${c%%:*}"
    COORD_RES="$COORD_RES,$r"
  done
  COORD_RES="${COORD_RES#,}"
  MN_POCK="(:$COORD_RES)|@MN"
fi
# 孤立 SASA 剥除掩码（配偶体不存在时不带前缀逗号）
STRIP_PROT="${PROT_P},:WAT,Na+,Cl-"; STRIP_PROT="${STRIP_PROT#,}"
STRIP_DNA="${DNA_P},:WAT,Na+,Cl-";   STRIP_DNA="${STRIP_DNA#,}"

# ---------------- an0：基础量（heredoc 条件拼装，输出 RES 顶层） ----------------
banner "an0_basics: RMSD/Rg/RMSF/SASA/PCA/Mn 配位时序"
AN0="$RAW/an0_basics.in"
{
  echo "parm $TOP_DIR/v.parm7"
  printf '%s\n' "$TRAJIN_BLOCK"
  echo "autoimage"
  [ "$PROT_HAS" = "1" ] && [ "$DNA_HAS" = "1" ] && echo "autoimage anchor $PROT_P"
  if [ "$PROT_HAS" = "1" ]; then
    cat <<EOF
rms PROT_BB out $RES_DIR/rmsd_prot.dat ${PROT_P}@N,CA,C,O first
rmsf PROT out $RES_DIR/rmsf_prot.dat ${PROT_P}@N,CA,C,O byres
radgyr PROT_RG ${PROT_P} out $RES_DIR/rg_prot.dat
surf ${PROT_P} out $RES_DIR/sasa_prot.dat
EOF
  fi
  if [ "$DNA_HAS" = "1" ]; then
    cat <<EOF
rms DNA_BB out $RES_DIR/rmsd_dna.dat ${DNA_P}@P,O5',C5',C4',O4',C3',O3' first
rmsf DNA out $RES_DIR/rmsf_dna.dat ${DNA_P} byres
radgyr DNA_RG ${DNA_P} out $RES_DIR/rg_dna.dat
surf ${DNA_P} out $RES_DIR/sasa_dna.dat
EOF
  fi
  if [ "$PROT_HAS" = "1" ] && [ "$DNA_HAS" = "1" ]; then
    echo "surf ${PROT_P},${DNA_P} out $RES_DIR/sasa_cplx.dat"
  elif [ "$PROT_HAS" = "1" ]; then
    echo "surf ${PROT_P} out $RES_DIR/sasa_cplx.dat"
  else
    echo "surf ${DNA_P} out $RES_DIR/sasa_cplx.dat"
  fi
  echo "run"
  # 孤立 SASA（strip+parmstrip 同步剥除；必须分别载入原始拓扑，cpptraj 不可回滚）
  if [ "$PROT_HAS" = "1" ]; then
    cat <<EOF
parm $TOP_DIR/v.parm7
trajin ${TRJ1} 1 last 1
autoimage
strip $STRIP_DNA
parmstrip $STRIP_DNA
surf ${PROT_P} out $RES_DIR/sasa_prot_iso.dat
run
EOF
  fi
  if [ "$DNA_HAS" = "1" ]; then
    cat <<EOF
parm $TOP_DIR/v.parm7
trajin ${TRJ1} 1 last 1
autoimage
strip $STRIP_PROT
parmstrip $STRIP_PROT
surf ${DNA_P} out $RES_DIR/sasa_dna_iso3.dat
run
EOF
  fi
  # PCA（蛋白 CA 全模式投影：pca_all=全部模式, pca_proj=前 3）
  if [ "$PROT_HAS" = "1" ]; then
    NCA_HI="${PROT_RANGE#*-}"; NCA_LO="${PROT_RANGE%-*}"
    NCA=$((NCA_HI - NCA_LO + 1)); N_MODE=$((NCA * 3))
    cat <<EOF
parm $TOP_DIR/v.parm7
trajin ${TRJ1} 1 last 1
rms FIT ${PROT_P}@CA first
matrix covar name PCM ${PROT_P}@CA
diagmatrix PCM out $RAW/pca_evecs.dat vecs 0
projection PALL $RAW/pca_evecs.dat ${PROT_P}@CA out $RES_DIR/pca_all.dat beg 1 end $N_MODE
projection P3   $RAW/pca_evecs.dat ${PROT_P}@CA out $RES_DIR/pca_proj.dat beg 1 end 3
run
EOF
  fi
  # Mn 配位距离时序（输出固定名 mn_H71/H75/E102a/E102b —— plot1 fig08 约定，见 README）
  if [ "$METAL" = "yes" ] && [ "${MN_RES:-}" != "-" ] && [ "${COORD1:-}" != "-" ]; then
    cat <<EOF
parm $TOP_DIR/v.parm7
trajin ${TRJ1} 1 last 1
EOF
    i=0
    for nm in H71 H75 E102a E102b; do
      i=$((i+1))
      c=$(eval "echo \${COORD$i:-}")
      [ "$c" = "-" ] && continue
      echo "distance mn_$nm :${MN_RES}@MN :${c%%:*}@${c##*:} out $RES_DIR/mn_$nm.dat"
    done
    echo "run"
  fi
} > "$AN0"
exec_cmd "'$CPPTRAJ' -i '$AN0' > '$RAW/an0_basics.log' 2>&1" || { tail -25 "$RAW/an0_basics.log" >&2; die "an0 采集失败"; }

# ---------------- an1..an5：专项量（模板渲染到 RAW） ----------------
for an in an1_struct an2b_contacts an2c_pock an2d_ee an4_dynamics an4b_dccm_fit an5_rdf; do
  # 依赖检查
  case "$an" in
    an1_struct|an4b_dccm_fit) [ "$PROT_HAS" = "1" ] || { info "跳过 $an（无蛋白）"; continue; } ;;
    an2b_contacts) [ "$PROT_HAS" = "1" ] && [ "$DNA_HAS" = "1" ] || { info "跳过 $an（需蛋白+DNA 界面）"; continue; } ;;
    an2d_ee|an4_dynamics) [ "$DNA_HAS" = "1" ] || { info "跳过 $an（无 DNA）"; continue; } ;;
    an2c_pock|an5_rdf) [ "$METAL" = "yes" ] && [ -n "${MN_POCK:-}" ] || { info "跳过 $an（无金属或无配位原子）"; continue; } ;;
  esac
  # an2c/an5 需要 Mn；an4/2d 需要 DNA 起止
  DNA_LO=""; DNA_HI=""; DNA_MID=""
  if [ "$DNA_HAS" = "1" ]; then
    DNA_LO="${DNA_RANGE%-*}"; DNA_HI="${DNA_RANGE#*-}"
    DNA_MID=$((DNA_LO + (DNA_HI - DNA_LO + 1) / 2))
  fi
  render "$CPPTRAJ_DIR/$an.in" "$RAW/$an.in" \
    PARM7="$TOP_DIR/v.parm7" TRAJIN_BLOCK="$TRAJIN_BLOCK" \
    PROT_RANGE="$PROT_RANGE" DNA_RANGE="$DNA_RANGE" \
    DNA_LO="$DNA_LO" DNA_MID="$DNA_MID" DNA_HI="$DNA_HI" \
    MN_RES="${MN_RES:--}" MN_POCK="${MN_POCK:--}" POCK_MASK="${MN_POCK:--}" \
    RAW="$RAW" RES="$RES_DIR"
  info "cpptraj $an"
  exec_cmd "'$CPPTRAJ' -i '$RAW/$an.in' > '$RAW/$an.log' 2>&1" || { tail -20 "$RAW/$an.log" >&2; die "$an 失败"; }
done

# an3_hbdist：体系特异关键氢键（config HBOND_PAIRS 定义），空则跳过
if [ -n "${HBOND_PAIRS:-}" ]; then
  HB_LINES=""
  IFS='|' read -ra pairs <<< "$HBOND_PAIRS"
  for pr in "${pairs[@]}"; do
    name="${pr%%=*}"; rest="${pr#*=}"
    a1="${rest%%=*}"; a2="${rest#*=}"
    # a1 = "146:OG" a2 = "271:O6"
    HB_LINES="$HB_LINES
distance hb_$name :${a1%%:*}@${a1##*:} :${a2%%:*}@${a2##*:} out $RAW/hbd_$name.dat"
  done
  render "$CPPTRAJ_DIR/an3_hbdist.in" "$RAW/an3_hbdist.in" \
    PARM7="$TOP_DIR/v.parm7" TRAJIN_BLOCK="$TRAJIN_BLOCK" HB_LINES="$HB_LINES"
  info "cpptraj an3_hbdist（${#pairs[@]} 条键）"
  exec_cmd "'$CPPTRAJ' -i '$RAW/an3_hbdist.in' > '$RAW/an3_hbdist.log' 2>&1" || { tail -20 "$RAW/an3_hbdist.log" >&2; die "an3 失败"; }
else
  info "HBOND_PAIRS 为空：跳过 an3（fig22 需在 config 定义关键键）"
fi

# 清理中间件
rm -f "$RAW/pca_evecs.dat"
ok "采集完成: $RES_DIR（顶层 .dat）+ cpptraj_raw/"
