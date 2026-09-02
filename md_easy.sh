#!/bin/bash
# ============================================================================
# md_easy.sh — md-amber-figure-pipeline 傻瓜式总入口（一条命令跑完整条 Amber MD）
#
#   bash md_easy.sh                 交互问答（全部回车=默认，也能跑通无金属蛋白体系）
#   bash md_easy.sh -c my.env       非交互：用已有配置跑全流程（正式模拟推荐）
#   bash md_easy.sh -c my.env -f    强制重跑（忽略断点）
#   bash md_easy.sh -c my.env --from 3    从阶段 3 续跑到结束
#   bash md_easy.sh -c my.env -s 3       只跑阶段 3
#   bash md_easy.sh -c my.env --dry-run  只做环境体检+预览命令，不真正执行
#   bash md_easy.sh -c my.env --env-only 只问答/校验配置后退出
#
# 阶段: 0 环境体检 → 1 tleap 拓扑 → 2 Mn/DISANG 约束 → 3 MD(min1→min2→heat→equil→prod)
#      → 4 cpptraj 数据采集 → 5 MM-PBSA(可选) → 6 一键 30 图
# 产物: $PWD/topol/$SYSTEM  $PWD/out/$SYSTEM  $PWD/results/$SYSTEM
# 断点: 每阶段成功写 .md_easy/$SYSTEM/stageN.done，重跑自动跳过已完成阶段
# ============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/pipeline/common.sh"
N_STAGES=6

usage() {
  sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'
  exit "${1:-0}"
}

# ---------- 解析命令行 ----------
CONFIG=""; RUN_STAGES=""; FROM=""; FORCE=0; DRYRUN=0; ENV_ONLY=0
while [ $# -gt 0 ]; do
  case "$1" in
    -c|--config) CONFIG="${2:-}"; shift 2 ;;
    -s|--stage)  RUN_STAGES="${2:-}"; shift 2 ;;
    --from)      FROM="${2:-}"; shift 2 ;;
    -f|--force)  FORCE=1; shift ;;
    --dry-run)   DRYRUN=1; shift ;;
    --env-only)  ENV_ONLY=1; shift ;;
    -h|--help)   usage 0 ;;
    *) die "未知参数: $1 （-h 看帮助）" ;;
  esac
done

# ---------- 配置来源 ----------
gen_env_interactive() { # 交互问答生成 project.env（默认值=回车可直通）
  banner "交互配置（回车接受 [默认值]）"
  printf '  输入复合物 PDB 路径: '; read -r p
  [ -n "$p" ] || die "PDB 路径必填（可拖文件进终端）"
  p="$(echo "$p" | sed 's/^["'"'"']//; s/["'"'"']$//')"
  [ -f "$p" ] || die "PDB 不存在: $p"
  local def_sys; def_sys="$(basename "$p" .pdb)"
  printf '  体系名 [%s]: ' "$def_sys"; read -r s; s="${s:-$def_sys}"
  printf '  金属离子 (Mn=含Mn2+体系 / none / auto=自动探测) [auto]: '; read -r m; m="${m:-auto}"
  printf '  生产时长 ns [150]: '; read -r n; n="${n:-150}"
  printf '  重复数 [1]: '; read -r r; r="${r:-1}"
  printf '  跑 MM-PBSA 逐残基分解 (y/N): '; read -r mm; mm="${mm:-n}"
  printf '  一键出 30 图 (Y/n): '; read -r fg; fg="${fg:-y}"
  case "$mm" in y|Y|yes|YES) mm=yes;; *) mm=no;; esac
  case "$fg" in n|N|no|NO) fg=no;; *) fg=yes;; esac
  cat > "$RUN_DIR/project.env" <<EOF
# md-amber-figure-pipeline 项目配置（由 md_easy.sh 生成）
SYSTEM="$s"
INPUT_PDB="$p"
METAL="$m"
PROD_NS="$n"
NREP="$r"
DO_MMPBSA="$mm"
DO_FIGS="$fg"
EOF
  ok "已生成 $RUN_DIR/project.env"
}

RUN_DIR="$PWD"

if [ -n "$CONFIG" ]; then
  [ -f "$CONFIG" ] || die "配置不存在: $CONFIG"
  RUN_DIR="$(cd "$(dirname "$CONFIG")" && pwd)"
  source "$CONFIG"
elif [ -f "$RUN_DIR/project.env" ] && [ "$FORCE" = "0" ] && [ -z "$FROM" ] && [ -z "$RUN_STAGES" ]; then
  # 普通全跑且已有配置 → 沿用（-f / --from / -s 也沿用但需走下方强制分支保证不丢）
  :
elif [ "$FORCE" = "1" ] || [ -n "$FROM" ] || [ -n "$RUN_STAGES" ] || [ -n "$CONFIG" ]; then
  [ -f "$RUN_DIR/project.env" ] && source "$RUN_DIR/project.env"
else
  # 无配置 → 交互生成（沿用 or 重问答）
  if [ -f "$RUN_DIR/project.env" ]; then
    printf '已存在 project.env，沿用？(y=沿用 / n=重新问答) [y]: '; read -r u
    if [ "$u" = "n" ] || [ "$u" = "N" ]; then gen_env_interactive; fi
    source "$RUN_DIR/project.env"
  else
    gen_env_interactive
    source "$RUN_DIR/project.env"
  fi
fi
export MDEASY_RUN_DIR="$RUN_DIR"

# ---------- 阶段解析 ----------
STAGE_LIST="0 1 2 3 4 5 6"
if [ -n "$FROM" ]; then
  case "$FROM" in [0-6]) ;; *) die "--from 取值 0-6" ;; esac
  STAGE_LIST="$(seq "$FROM" 6 | tr '\n' ' ')"
elif [ -n "$RUN_STAGES" ]; then
  STAGE_LIST=""
  IFS=',' read -ra parts <<< "$RUN_STAGES"
  for p in "${parts[@]}"; do
    case "$p" in
      [0-6]) STAGE_LIST="$STAGE_LIST $p" ;;
      [0-6]-[0-6]) lo="${p%-*}"; hi="${p#*-}"; STAGE_LIST="$STAGE_LIST $(seq "$lo" "$hi")" ;;
      *) die "非法阶段: $p（取 0-6，可用逗号/区间，如 1,3-5）" ;;
    esac
  done
  STAGE_LIST="$(echo "$STAGE_LIST" | tr ' ' '\n' | sort -n -u | tr '\n' ' ')"
fi
export MDEASY_DRYRUN="$DRYRUN"

# ---------- 初始化（探测环境但暂不 die：stage0 统一体检更友好）----------
load_defaults
[ -n "${INPUT_PDB:-}" ] && [ -f "$INPUT_PDB" ] && init_workdirs
if [ "$ENV_ONLY" = "1" ]; then
  banner "配置预览"
  env | grep -E "^(SYSTEM|INPUT_PDB|PROT_FF|DNA_FF|WAT_FF|BOX_DIST|METAL|FRCMOD|PROD_NS|NREP|TEMP|DO_MMPBSA|DO_FIGS|HBOND_PAIRS)=" | sort
  exit 0
fi

# ---------- 工作目录就绪（-s/--from 无 config 场景同样保证变量） ----------
[ -z "${MDEASY_RUN_DIR:-}" ] && MDEASY_RUN_DIR="$PWD"
init_workdirs
cd "$RUN_DIR"
if [ -z "${INPUT_PDB:-}" ] || [ ! -f "$INPUT_PDB" ]; then
  die "INPUT_PDB 缺失或不存在（当前目录: $PWD）: ${INPUT_PDB:-<未设置>} —— 请用 -c 指定配置或重新问答"
fi
[ "$MDEASY_DRYRUN" = "1" ] && info "DRY-RUN 模式：仅体检与预览命令"

# ---------- 执行 ----------
for s in $STAGE_LIST; do
  case "$s" in
    0) f="$HERE/pipeline/stage0_env.sh"; t="env";;
    1) f="$HERE/pipeline/stage1_tleap.sh"; t="tleap";;
    2) f="$HERE/pipeline/stage2_restraints.sh"; t="restraints";;
    3) f="$HERE/pipeline/stage3_md.sh"; t="md";;
    4) f="$HERE/pipeline/stage4_cpptraj.sh"; t="cpptraj";;
    5) f="$HERE/pipeline/stage5_mmpbsa.sh"; t="mmpbsa";;
    6) f="$HERE/pipeline/stage6_figs.sh"; t="figs";;
  esac
  if [ "$FORCE" = "0" ] && stage_mark "$s" "$t"; then
    ok "阶段 $s ($t) 已完成，跳过（-f 强制重跑）"
    continue
  fi
  if [ "$DRYRUN" = "1" ] && [ "$s" -gt 0 ]; then
    banner "阶段 $s ($t) [dry-run 预览]"
    bash "$f" || true
    continue
  fi
  banner "阶段 $s/$N_STAGES — $t"
  if bash "$f"; then
    stage_done "$s" "$t"
    ok "阶段 $s ($t) 完成"
  else
    err "阶段 $s ($t) 失败 —— 日志见 $LOG_DIR/ ；修复后重跑同命令即可断点续跑"
    exit 1
  fi
done

banner "全部完成"
say "  拓扑    : $TOP_DIR"
say "  轨迹    : $OUT_DIR"
say "  分析/图 : $RES_DIR"
say "重跑续跑  : bash $0 -c <config> --from 0 （或直接同命令）"
