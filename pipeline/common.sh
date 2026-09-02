# ============================================================================
# common.sh — md-amber-figure-pipeline 公共库（被 md_easy.sh 与 pipeline/*.sh source）
#   职责: 仓库定位 / 配色日志 / 配置默认值+加载 / 派生路径 / 引擎探测 / 命令执行
# ============================================================================
[ -n "${MD_COMMON_LOADED:-}" ] && return 0
MD_COMMON_LOADED=1

# ---- 仓库定位（允许被软链/拷贝；以本文件真实路径为准）----
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$PIPELINE_DIR")"
TOOLS_DIR="$REPO_ROOT/tools"
TPL_DIR="$REPO_ROOT/templates"
CPPTRAJ_DIR="$REPO_ROOT/cpptraj"
SCRIPTS_DIR="$REPO_ROOT/scripts"
FF_DIR="$REPO_ROOT/ff"
PYTHON="${PYTHON:-python3}"

# ---- 配色（仅 TTY）----
if [ -t 1 ]; then
  C_OK=$'\033[32m'; C_WARN=$'\033[33m'; C_ERR=$'\033[31m'; C_INF=$'\033[36m'; C_B=$'\033[1m'; C0=$'\033[0m'
else
  C_OK=""; C_WARN=""; C_ERR=""; C_INF=""; C_B=""; C0=""
fi
say()   { printf '%s\n' "$*"; }
ok()    { printf '%s[ OK ]%s %s\n' "$C_OK" "$C0" "$*"; }
warn()  { printf '%s[WARN]%s %s\n' "$C_WARN" "$C0" "$*" >&2; }
err()   { printf '%s[FAIL]%s %s\n' "$C_ERR" "$C0" "$*" >&2; }
die()   { err "$*"; exit 1; }
info()  { printf '%s[----]%s %s\n' "$C_INF" "$C0" "$*"; }
banner(){ printf '\n%s==== %s ====%s\n' "$C_B" "$*" "$C0"; }
need()  { command -v "$1" >/dev/null 2>&1 || die "缺少必需命令: $1"; }

# ---- 命令执行（dry-run 感知；接受整条命令行字符串，bash -c 执行）----
exec_cmd() {
  local cmd="$1"
  if [ "${MDEASY_DRYRUN:-0}" = "1" ]; then
    info "(dry-run) $cmd"
    return 0
  fi
  bash -c "$cmd" || die "命令失败 → $cmd"
}
# 不中断失败的执行（用于"尝试 GPU 失败退 CPU"）
exec_cmd_nofail() {
  local cmd="$1"
  if [ "${MDEASY_DRYRUN:-0}" = "1" ]; then
    info "(dry-run) $cmd"
    return 0
  fi
  bash -c "$cmd"
}

# ---- 模板渲染（render_tmpl.py 逐占位符替换，规避 sed 对 mask 中 | & 的坑）----
render() { # render <模板> <输出> [KEY=VAL ...]
  "$PYTHON" "$TOOLS_DIR/render_tmpl.py" "$@" || die "模板渲染失败: $1"
}

# ---- 阶段标记（断点续跑）----
stage_mark() { # stage_mark <N> <name>
  local dir="${MARK_DIR:?MARK_DIR 未设置}"
  mkdir -p "$dir"
  [ -f "$dir/stage${1}_${2}.done" ]
}
stage_done() {
  mkdir -p "${MARK_DIR:?}"
  : > "$MARK_DIR/stage${1}_${2}.done"
}
stage_clear_from() { # 清除 >=N 的所有阶段标记（--from 时用）
  local n="$1"
  [ -d "${MARK_DIR:-}" ] || return 0
  for f in "$MARK_DIR"/stage*.done; do
    [ -f "$f" ] || continue
    local tag; tag="$(basename "$f" .done)"; tag="${tag#stage}"
    local num="${tag%%_*}"
    if [ "$num" -ge "$n" ] 2>/dev/null; then rm -f "$f"; fi
  done
}

# ---- 默认配置（project.env 未设项全部兜底）----
load_defaults() {
  export SYSTEM="${SYSTEM:-sys1}"
  export INPUT_PDB="${INPUT_PDB:-}"
  export PROT_FF="${PROT_FF:-ff19SB}"
  export DNA_FF="${DNA_FF:-OL15}"
  export WAT_FF="${WAT_FF:-TIP3P}"
  export BOX_DIST="${BOX_DIST:-12.0}"
  export METAL="${METAL:-auto}"
  export FRCMOD="${FRCMOD:-}"
  export NREP="${NREP:-1}"
  export PROD_NS="${PROD_NS:-150}"
  export TEMP="${TEMP:-300.0}"
  export MIN_WT="${MIN_WT:-25.0}"
  export HEAT_WT="${HEAT_WT:-10.0}"
  export EQUIL_WT="${EQUIL_WT:-2.0}"
  export ENGINE="${ENGINE:-auto}"
  export DO_MMPBSA="${DO_MMPBSA:-no}"
  export DO_FIGS="${DO_FIGS:-yes}"
  export PLOT_DT_NS="${PLOT_DT_NS:-0.1}"
  export HBOND_PAIRS="${HBOND_PAIRS:-}"
}

# ---- 工作区布局（默认 RUN_DIR=调用目录；可用 MDEASY_RUN_DIR 覆盖）----
init_workdirs() {
  load_defaults
  RUN_DIR="${MDEASY_RUN_DIR:-$PWD}"
  LOG_DIR="$RUN_DIR/logs/$SYSTEM"
  TOP_DIR="$RUN_DIR/topol/$SYSTEM"
  OUT_DIR="$RUN_DIR/out/$SYSTEM"
  RES_DIR="$RUN_DIR/results/$SYSTEM"
  MARK_DIR="$RUN_DIR/.md_easy/$SYSTEM"
  mkdir -p "$LOG_DIR" "$TOP_DIR" "$OUT_DIR" "$RES_DIR" "$MARK_DIR"
}

# ---- 引擎探测（Amber / cpptraj / GPU）----
detect_amber() {
  if [ -z "${AMBERHOME:-}" ] || [ ! -x "$AMBERHOME/bin/tleap" ]; then
    local tp; tp="$(command -v tleap 2>/dev/null || true)"
    if [ -n "$tp" ]; then
      # .../bin/tleap -> AMBERHOME = 上上级目录（去掉 bin/tleap）
      export AMBERHOME="$(cd "$(dirname "$tp")/.." && pwd)"
    fi
  fi
  if [ -z "${AMBERHOME:-}" ] || [ ! -x "$AMBERHOME/bin/tleap" ]; then
    die "未找到 AmberTools（tleap）。请 source 你的 amber 环境后重试，或设置 AMBERHOME。"
  fi
  export PMEMD_CUDA="$AMBERHOME/bin/pmemd.cuda"
  export PMEMD_CPU="$AMBERHOME/bin/pmemd"
  export TLEAP="$AMBERHOME/bin/tleap"
  # cpptraj: 优先 AMBERHOME/bin，其次 PATH
  if [ -x "$AMBERHOME/bin/cpptraj" ]; then
    export CPPTRAJ="$AMBERHOME/bin/cpptraj"
  else
    CPPTRAJ="$(command -v cpptraj 2>/dev/null || true)"
    [ -n "$CPPTRAJ" ] || die "未找到 cpptraj（AMBERHOME/bin 与 PATH 都没有）"
    export CPPTRAJ
  fi
  # 引擎选择
  export MD_ENGINE="$PMEMD_CPU"
  if [ "$ENGINE" = "auto" ] || [ "$ENGINE" = "gpu" ]; then
    if [ -x "$PMEMD_CUDA" ]; then MD_ENGINE="$PMEMD_CUDA"; fi
  fi
  ok "AMBERHOME=$AMBERHOME"
  info "MD 引擎: $MD_ENGINE   (GPU 失败将自动退回 $PMEMD_CPU)"
}

# ---- 金属判定 ----
# PDB 含 MN HETATM → METAL 归一化为 yes/no；METAL=auto 由 PDB 决定
resolve_metal() {
  if [ "$METAL" = "auto" ]; then
    if grep -q "^HETATM.* MN " "$INPUT_PDB" 2>/dev/null || grep -qE "^HETATM.{8}MN " "$INPUT_PDB" 2>/dev/null; then
      export METAL="yes"
    else
      export METAL="no"
    fi
  fi
  if [ "$METAL" = "yes" ]; then
    export FRCMOD="${FRCMOD:-$FF_DIR/mn_cm12-6.frcmod}"
    [ -f "$FRCMOD" ] || die "金属 frcmod 不存在: $FRCMOD（METAL=yes 时必须有）"
  fi
  info "METAL=$METAL"
}

# ---- 阶段 banner 包装 ----
stage_header() { # stage_header <N> <名称>  （N_STAGES 由 md_easy.sh 定义; 单跑 stage 脚本时兜底 6）
  banner "阶段 $1/${N_STAGES:-6} — $2"
}
