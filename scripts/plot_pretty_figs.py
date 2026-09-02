#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PprI(WT)·S1_ssDNA·Mn  150 ns 显式溶剂 MD —— 知乎推文级 7 图套件 (+1 机制图)

数据来源: CHPC 4090 真实 MD (prod.nc, 1500 帧 × 100 ps = 150 ns)
  蛋白 1-254 (seq 编号 = PDB 作者号 - 21), DNA 255-278 (S1_G17, 24 nt), Mn
产出: results/WT__S1/figures_pretty/fig01..fig08

作者: WorkBuddy   日期: 2026-09-02
"""
import os, re, io
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FuncFormatter
from matplotlib import font_manager
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from scipy.stats import gaussian_kde
from scipy.ndimage import gaussian_filter

# ----------------------------------------------------------------- 路径
HERE = os.path.dirname(os.path.abspath(__file__))


def _norm_msys(p):
    """Git Bash(MSYS) 路径 /a/Data/... → A:/Data/...（Windows 原生 python 可读）"""
    m = re.match(r"^/([a-zA-Z])/(.*)$", p)
    return f"{m.group(1).upper()}:/{m.group(2)}" if m else p


# 傻瓜化: 允许 env 覆盖（md_easy.sh stage6 传入 MDEASY_RES=<abs results/$SYSTEM>）
_RES_OVERRIDE = _norm_msys(os.environ.get("MDEASY_RES") or
                           os.path.join(HERE, "results", "WT__S1"))
CASE = os.environ.get("MDEASY_CASE") or "WT__S1"
RES = _RES_OVERRIDE
OUT = os.path.join(RES, "figures_pretty")
os.makedirs(OUT, exist_ok=True)
# env 覆盖模式: 非 PprI 体系时关闭 PprI 残基注释（fig02 等仍出图但省注释）
ENV_OVERRIDE = bool(os.environ.get("MDEASY_RES"))
if ENV_OVERRIDE:
    print(f"[通用模式] MDEASY_RES={RES}  关闭 PprI 特异残基注释")


def _ttl(desc):
    """通用 suptitle 前缀: 非 PprI 体系(MDEASY_RES)时用 CASE, 原版保留 PprI 标签"""
    dur = f"{NFR * DT_NS:.0f} ns" if "NFR" in globals() else "MD"
    if ENV_OVERRIDE:
        return f"{CASE}  {dur} MD — {desc}"
    return f"PprI(WT)·S1-ssDNA·Mn$^{{2+}}$   {dur} MD — {desc}"


# ----------------------------------------------------------------- 全局样式
YAHEI = "Microsoft YaHei"
for f in font_manager.fontManager.ttflist:
    if f.name == YAHEI:
        break

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": [YAHEI, "DejaVu Sans"],
    "axes.unicode_minus": False,
    "axes.linewidth": 1.3,
    "axes.edgecolor": "#2B2B2B",
    "axes.labelcolor": "#1A1A1A",
    "axes.labelsize": 12.5,
    "axes.titlesize": 14.5,
    "axes.titleweight": "bold",
    "axes.titlepad": 10,
    "xtick.labelsize": 10.5,
    "ytick.labelsize": 10.5,
    "xtick.major.size": 4.5,
    "ytick.major.size": 4.5,
    "xtick.major.width": 1.3,
    "ytick.major.width": 1.3,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.color": "#2B2B2B",
    "ytick.color": "#2B2B2B",
    "legend.frameon": False,
    "legend.fontsize": 10.5,
    "legend.handlelength": 1.8,
    "legend.columnspacing": 1.2,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# ----------------------------------------------------------------- 配色
C_PROT = "#2B6CB0"      # 蛋白 深蓝
C_DNA = "#E8623C"       # DNA 橙红
C_CPLX = "#3E8E6E"      # 复合物 绿
C_ACC = "#F0A202"       # 强调 琥珀
C_PUR = "#7B4EA8"       # 紫
C_TEAL = "#1B9AAA"      # 青
GREY = "#98A2AB"
INK = "#1A1A1A"
GRIDC = "#D8DEE3"

PALETTE = [C_PROT, C_DNA, C_CPLX, C_ACC, C_PUR, C_TEAL]

# ----------------------------------------------------------------- 残基注释 (seq 编号 = PDB - 21)
OFF = 21
KEY_RES = {
    64: ("R85", "Patch1"), 67: ("F88", "Core1"), 71: ("H92", "HEXXH"),
    75: ("H96", "HEXXH"), 102: ("E123", "HEXXH"), 106: ("N127", "Core1"),
    146: ("S167", "Core1"), 149: ("Y170", "Core1"),
    186: ("R207", "Patch1"), 196: ("Y217", "Core2"), 199: ("R220", "Patch3"),
    229: ("R250", "Patch3"), 230: ("S251", "Patch3"), 232: ("R253", "Core2"),
    234: ("M255", "Core2"), 246: ("R267", "Patch1"),
}
HELIX_CAT = (62, 80)      # 催化螺旋 (PDB 83-101, 含 HEXXH)
LOOP_MOD = (168, 181)     # 8SLN 无序环 (PDB 189-202), 建模补出

GROUP_COLOR = {"Core1": C_PROT, "Core2": C_PUR, "HEXXH": C_ACC,
               "Patch1": C_TEAL, "Patch3": C_DNA}


# ----------------------------------------------------------------- 工具
def read_dat(path, skip=1):
    """读 cpptraj .dat: 首行为 '#Frame ...' 表头"""
    return np.loadtxt(path, skiprows=skip)


def roll_mean(y, w=10):
    """居中滚动平均, 返回与原数组等长 (两端用边缘值填充)"""
    y = np.asarray(y, float)
    if len(y) < w:
        return y.copy()
    k = np.ones(w) / w
    ypad = np.concatenate([np.full(w // 2, y[0]), y, np.full(w - 1 - w // 2, y[-1])])
    return np.convolve(ypad, k, mode="valid")


def style_ax(ax, gridy=True, gridx=False):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.3)
    ax.spines["bottom"].set_linewidth(1.3)
    if gridy:
        ax.grid(axis="y", ls="--", lw=0.7, color=GRIDC, alpha=0.75, zorder=0)
    if gridx:
        ax.grid(axis="x", ls="--", lw=0.7, color=GRIDC, alpha=0.75, zorder=0)
    ax.set_axisbelow(True)


def panel_label(ax, txt, x=-0.02, y=1.06):
    if hasattr(ax, "text2D"):            # Axes3D: text() 需要 (x,y,z)
        ax.text2D(x, y, txt, transform=ax.transAxes, fontsize=17,
                  fontweight="bold", color=INK, va="bottom", ha="right")
        return
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=17, fontweight="bold",
            color=INK, va="bottom", ha="right")


def add_marginal(fig, gs_cell, data, color, xlabel=None, bins=42,
                 orient="right", share_ylim=None, lw=1.6):
    """在时间序主图右侧加 KDE + 直方图边际分布"""
    ax = fig.add_subplot(gs_cell)
    d = np.asarray(data, float)
    d = d[~np.isnan(d)]
    kde = gaussian_kde(d) if len(d) > 5 else None
    counts, edges = np.histogram(d, bins=bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    ax.barh(centers, counts, height=(edges[1] - edges[0]) * 0.92,
            color=color, alpha=0.22, edgecolor="none", zorder=1)
    lo, hi = (share_ylim if share_ylim else (d.min(), d.max()))
    pad = 0.12 * (hi - lo + 1e-9)
    ys = np.linspace(lo - pad, hi + pad, 300)
    if kde is not None:
        ax.plot(kde(ys), ys, color=color, lw=lw, zorder=3)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_xlim(0, max(kde(ys).max() if kde is not None else 1,
                       counts.max()) * 1.22)
    ax.axis("off")
    if xlabel:
        ax.text(0.5, 1.012, xlabel, transform=ax.transAxes, ha="center",
                va="bottom", fontsize=9.5, color="#5A6670")
    return ax


def save(fig, name):
    p = os.path.join(OUT, name)
    fig.savefig(p, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  ->", name)


def stat_box(ax, txt, loc="upper left", fs=9.2):
    ax.text(0.015 if "left" in loc else 0.985, 0.975, txt,
            transform=ax.transAxes, ha="left" if "left" in loc else "right",
            va="top", fontsize=fs, color="#3C4852", linespacing=1.55,
            bbox=dict(boxstyle="round,pad=0.45", fc="white", ec=GRIDC, lw=1.0, alpha=0.94))


# ================================================================= 载入数据
print("[load] reading cpptraj data ...")
# 必需核心数据快速检查 —— 缺哪个直接提示, 避免深埋在 read_dat 里的 FileNotFoundError
for _req in ("rmsd_prot.dat", "rmsf_prot.dat", "rg_prot.dat", "sasa_prot.dat",
             "sasa_prot_iso.dat", "sasa_cplx.dat", "pca_proj.dat", "pca_all.dat"):
    if not os.path.isfile(os.path.join(RES, _req)):
        raise SystemExit(f"[error] 缺必需数据 {os.path.join(RES, _req)}\n"
                         f"  → 请先运行 stage4 (bash md_easy.sh -c <config> --from 4)\n"
                         f"  → 或检查 MDEASY_RES 是否指向含 cpptraj 产物的 results/<SYSTEM>/")
D = {}
D["rmsd_prot"] = read_dat(os.path.join(RES, "rmsd_prot.dat"))[:, 1]
D["rmsd_dna"] = read_dat(os.path.join(RES, "rmsd_dna.dat"))[:, 1]
rg_p = read_dat(os.path.join(RES, "rg_prot.dat"))
D["rg_prot"], D["rg_prot_max"] = rg_p[:, 1], rg_p[:, 2]
rg_d = read_dat(os.path.join(RES, "rg_dna.dat"))
D["rg_dna"], D["rg_dna_max"] = rg_d[:, 1], rg_d[:, 2]
D["rmsf_prot"] = read_dat(os.path.join(RES, "rmsf_prot.dat"))[:, 1]
D["rmsf_dna"] = read_dat(os.path.join(RES, "rmsf_dna.dat"))[:, 1]
D["sasa_prot"] = read_dat(os.path.join(RES, "sasa_prot.dat"))[:, 1]
D["sasa_dna"] = read_dat(os.path.join(RES, "sasa_dna.dat"))[:, 1]
D["sasa_cplx"] = read_dat(os.path.join(RES, "sasa_cplx.dat"))[:, 1]
D["sasa_prot_iso"] = read_dat(os.path.join(RES, "sasa_prot_iso.dat"))[:, 1]
D["sasa_dna_iso"] = read_dat(os.path.join(RES, "sasa_dna_iso3.dat"))[:, 1]
# Mn²⁺ 配位时序 —— 无金属体系文件不存在, 条件装载(fig08 依 D 键存在性出占位图)
for _mk, _mf in (("mn_h71", "mn_H71.dat"), ("mn_h75", "mn_H75.dat"),
                 ("mn_e102a", "mn_E102a.dat"), ("mn_e102b", "mn_E102b.dat")):
    _p = os.path.join(RES, _mf)
    if os.path.isfile(_p):
        D[_mk] = read_dat(_p)[:, 1]
if not any(k in D for k in ("mn_h71", "mn_h75", "mn_e102a", "mn_e102b")):
    print("[load] 无 Mn 配位数据 (metal-free) — fig08 将输出占位图")
pca = read_dat(os.path.join(RES, "pca_proj.dat"))     # 1500 x (frame, m1, m2, m3)
D["pc1"], D["pc2"], D["pc3"] = pca[:, 1], pca[:, 2], pca[:, 3]
print("[load] pca_all (762 modes) ...")
pca_all = np.loadtxt(os.path.join(RES, "pca_all.dat"), skiprows=1)
PC_ALL = pca_all[:, 1:]
eig = PC_ALL.var(axis=0, ddof=1)
eig_sorted = np.sort(eig)[::-1]
var_pct = 100.0 * eig_sorted / eig_sorted.sum()

NFR = len(D["rmsd_prot"])
DT_NS = 0.1                      # 100 ps / frame
T = np.arange(NFR) * DT_NS
print(f"[load] {NFR} frames, {T[-1]:.1f} ns ; "
      f"PC1 {var_pct[0]:.1f}% PC2 {var_pct[1]:.1f}% PC3 {var_pct[2]:.1f}%")

# DNA 残基编号 → 序列位置 (255..278 -> 1..24)
DNA_POS = np.arange(255, 255 + len(D["rmsf_dna"]))


# ================================================================= FIG 1  RMSD
print("[fig01] RMSD")
fig = plt.figure(figsize=(9.4, 4.5))
gs = fig.add_gridspec(1, 3, width_ratios=[5.4, 0.85, 0.85], wspace=0.06,
                      left=0.075, right=0.975, top=0.86, bottom=0.15)
ax = fig.add_subplot(gs[0, 0])
for key, col, lab in [("rmsd_prot", C_PROT, "Protein backbone"),
                      ("rmsd_dna", C_DNA, "ssDNA backbone")]:
    y = D[key]
    ax.plot(T, y, lw=0.85, color=col, alpha=0.32, zorder=2)
    ax.plot(T, roll_mean(y, 10), lw=2.0, color=col, label=lab, zorder=3)
style_ax(ax)
ax.set_xlabel("Time (ns)", fontsize=12.5)
ax.set_ylabel(r"C$\alpha$/backbone RMSD ($\rm \AA$)", fontsize=12.5)
ax.set_xlim(0, T[-1])
ax.xaxis.set_major_locator(MultipleLocator(25))
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.legend(loc="upper left", ncol=2, frameon=False, bbox_to_anchor=(0.005, 1.005))
for key, col in [("rmsd_prot", C_PROT), ("rmsd_dna", C_DNA)]:
    ax.axhline(np.nanmean(D[key]), ls=":", lw=1.2, color=col, alpha=0.65, zorder=1)
txt = (f"Protein  {np.nanmean(D['rmsd_prot']):.2f} ± {np.nanstd(D['rmsd_prot']):.2f} $\\rm\\AA$\n"
       f"ssDNA      {np.nanmean(D['rmsd_dna']):.2f} ± {np.nanstd(D['rmsd_dna']):.2f} $\\rm\\AA$")
stat_box(ax, txt, loc="upper right")
add_marginal(fig, gs[0, 1], D["rmsd_prot"], C_PROT, xlabel="Protein",
             share_ylim=(0, max(np.nanmax(D["rmsd_prot"]), np.nanmax(D["rmsd_dna"])) * 1.02))
add_marginal(fig, gs[0, 2], D["rmsd_dna"], C_DNA, xlabel="ssDNA",
             share_ylim=(0, max(np.nanmax(D["rmsd_prot"]), np.nanmax(D["rmsd_dna"])) * 1.02))
fig.suptitle(_ttl("构象稳定性 (RMSD)"),
             fontsize=14.5, fontweight="bold", y=0.975)
save(fig, "fig01_rmsd.png")


# ================================================================= FIG 2  RMSF
print("[fig02] RMSF")
fig = plt.figure(figsize=(12.2, 6.1))
gs = fig.add_gridspec(2, 2, width_ratios=[5.2, 1.0], height_ratios=[4.6, 0.62],
                      wspace=0.06, hspace=0.09,
                      left=0.062, right=0.978, top=0.885, bottom=0.125)
ax = fig.add_subplot(gs[0, 0])
res = np.arange(1, len(D["rmsf_prot"]) + 1)
ax.fill_between(res, 0, D["rmsf_prot"], color=C_PROT, alpha=0.16, zorder=1)
ax.plot(res, D["rmsf_prot"], lw=1.25, color=C_PROT, zorder=3, label="Protein C$\\alpha$ RMSF")
ax.plot(res, roll_mean(D["rmsf_prot"], 7), lw=2.0, color="#0F4C81", zorder=4)
style_ax(ax)
ax.set_xlim(1, len(res))
ax.set_ylim(0, max(np.nanmax(D["rmsf_prot"]) * 1.28, 6))
ax.set_ylabel(r"RMSF ($\rm \AA$)", fontsize=12.5)
ax.xaxis.set_major_locator(MultipleLocator(25))
ax.tick_params(labelbottom=False)
ax.legend(loc="upper right", frameon=False)

# 结构域/功能元件条带 (PprI 特异; env 覆盖模式下跳过 280-302 注释段)
axb = fig.add_subplot(gs[1, 0], sharex=ax)
axb.axis("off")
if ENV_OVERRIDE:
    axb.plot([1, len(res)], [0.30, 0.30], lw=1.2, color=INK,
             transform=axb.get_xaxis_transform(), clip_on=False)
    for r in range(1, len(res) + 1, 25):
        axb.text(r, 0.10, str(r), ha="center", va="top", fontsize=8, color="#6B7681",
                 transform=axb.get_xaxis_transform(), clip_on=False)
    axb.set_ylim(0, 1)
    axb.set_xticks([])
else:
  bands = [(HELIX_CAT[0], HELIX_CAT[1], C_ACC, "catalytic helix 83–101 (HEXXH)"),
           (LOOP_MOD[0], LOOP_MOD[1], GREY, "disordered loop 189–202 (modelled)")]
  for a, b, c, lab in bands:
      axb.add_patch(plt.Rectangle((a, 0.10), b - a, 0.34, color=c, alpha=0.85,
                                  transform=axb.get_xaxis_transform(), clip_on=False))
      axb.text((a + b) / 2, 0.52, lab, ha="center", va="bottom", fontsize=8.4, color="#4A5560")
  axb.set_ylim(0, 1)
  axb.plot([1, len(res)], [0.06, 0.06], lw=1.0, color=INK,
           transform=axb.get_xaxis_transform(), clip_on=False)
  for r in [1, 62, 80, 168, 181, 254]:
      axb.plot([r, r], [0.02, 0.10], lw=1.0, color=INK)
      axb.text(r, -0.08, str(r + OFF), ha="center", va="top", fontsize=8, color="#6B7681")

  # 关键残基标注 (精选子集 + 双层交错防重叠)
  ANNOT = [64, 67, 71, 75, 102, 106, 146, 149, 186, 196, 232, 246]
  for i, r in enumerate(sorted(ANNOT)):
      nm, grp = KEY_RES[r]
      if r > len(D["rmsf_prot"]):
          continue
      y0 = D["rmsf_prot"][r - 1]
      up = bool(i % 2)
      yy = y0 + (0.45 if up else 1.5)
      ax.plot([r], [y0], "o", ms=6.2, color=GROUP_COLOR[grp], mec="white", mew=1.3, zorder=6)
      ax.annotate(f"{nm}", xy=(r, y0), xytext=(r, yy), ha="center", va="bottom",
                  fontsize=8.6, fontweight="bold", color=GROUP_COLOR[grp], zorder=6,
                  arrowprops=dict(arrowstyle="-", lw=0.9, color=GROUP_COLOR[grp], alpha=0.85))
  axb.text(0.5, -0.62, "Residue (PDB author numbering)", transform=axb.transAxes,
           ha="center", va="top", fontsize=12.5, color=INK)
  axb.set_xticks([])

# 右侧: DNA 逐核苷酸 RMSF
axd = fig.add_subplot(gs[:, 1])
pos = np.arange(1, len(D["rmsf_dna"]) + 1)
cols = [C_DNA if D["rmsf_dna"][i] >= np.median(D["rmsf_dna"]) else "#F5B48A"
        for i in range(len(pos))]
axd.barh(pos, D["rmsf_dna"], color=cols, edgecolor="white", lw=0.6, height=0.78)
style_ax(axd, gridy=False, gridx=True)
axd.set_yticks([1, 5, 10, 15, 20, 24])
axd.set_ylim(0.2, len(pos) + 0.8)
axd.set_xlabel(r"RMSF ($\rm \AA$)", fontsize=11.5)
axd.set_ylabel("ssDNA position", fontsize=11.5)
axd.invert_yaxis()
axd.set_title("ssDNA (24 nt)", fontsize=12, fontweight="bold", pad=8)
handles = [plt.Line2D([], [], marker="o", ls="", color=c, label=g)
           for g, c in GROUP_COLOR.items()]
ax.legend(handles=handles, loc="upper left", fontsize=8.8, frameon=False,
          title="functional group", title_fontsize=8.8,
          bbox_to_anchor=(0.155, 1.005))
panel_label(ax, "A")
panel_label(axd, "B", x=-0.06, y=1.06)
fig.suptitle(_ttl("逐残基柔性 (RMSF)"),
             fontsize=14.5, fontweight="bold", y=0.985)
save(fig, "fig02_rmsf.png")


# ================================================================= FIG 3  Rg
print("[fig03] Rg")
fig = plt.figure(figsize=(9.4, 4.6))
gs = fig.add_gridspec(1, 2, width_ratios=[5.6, 1.0], wspace=0.06,
                      left=0.075, right=0.975, top=0.86, bottom=0.15)
ax = fig.add_subplot(gs[0, 0])
for key, col, lab in [("rg_prot", C_PROT, "Protein"), ("rg_dna", C_DNA, "ssDNA")]:
    y = D[key]
    ax.plot(T, y, lw=0.85, color=col, alpha=0.30, zorder=2)
    ax.plot(T, roll_mean(y, 10), lw=2.0, color=col, label=lab, zorder=3)
    ax.axhline(np.nanmean(y), ls=":", lw=1.2, color=col, alpha=0.6, zorder=1)
style_ax(ax)
ax.set_xlabel("Time (ns)", fontsize=12.5)
ax.set_ylabel(r"Radius of gyration, $R_g$ ($\rm \AA$)", fontsize=12.5)
ax.set_xlim(0, T[-1])
ax.xaxis.set_major_locator(MultipleLocator(25))
ax.legend(loc="upper right", frameon=False, ncol=2)
txt = (f"Protein  {np.nanmean(D['rg_prot']):.2f} ± {np.nanstd(D['rg_prot']):.2f} $\\rm\\AA$\n"
       f"ssDNA      {np.nanmean(D['rg_dna']):.2f} ± {np.nanstd(D['rg_dna']):.2f} $\\rm\\AA$")
stat_box(ax, txt, loc="upper left")
# 边际: 两个 KDE 叠加
axm = fig.add_subplot(gs[0, 1])
allv = np.concatenate([D["rg_prot"], D["rg_dna"]])
lo, hi = allv.min(), allv.max()
pad = 0.1 * (hi - lo)
ys = np.linspace(lo - pad, hi + pad, 300)
for key, col, lab in [("rg_prot", C_PROT, "Protein"), ("rg_dna", C_DNA, "ssDNA")]:
    d = D[key]
    kde = gaussian_kde(d)
    axm.plot(kde(ys), ys, color=col, lw=1.7, label=lab)
    cnt, edg = np.histogram(d, bins=40, density=True)
    axm.barh(0.5 * (edg[:-1] + edg[1:]), cnt, height=(edg[1] - edg[0]) * 0.9,
             color=col, alpha=0.18, edgecolor="none")
axm.set_ylim(lo - pad, hi + pad)
axm.axis("off")
axm.text(0.5, 1.012, "Density", transform=axm.transAxes, ha="center",
         va="bottom", fontsize=9.5, color="#5A6670")
fig.suptitle(_ttl("整体紧致度 (Rg)"),
             fontsize=14.5, fontweight="bold", y=0.975)
save(fig, "fig03_rg.png")


# ================================================================= FIG 4  SASA
print("[fig04] SASA")
fig = plt.figure(figsize=(10.4, 4.9))
gs = fig.add_gridspec(1, 2, width_ratios=[5.6, 1.0], wspace=0.06,
                      left=0.072, right=0.975, top=0.86, bottom=0.145)
ax = fig.add_subplot(gs[0, 0])
for key, col, lab in [("sasa_cplx", C_CPLX, "Complex"),
                      ("sasa_prot", C_PROT, "Protein"),
                      ("sasa_dna", C_DNA, "ssDNA")]:
    y = D[key]
    ax.plot(T, y, lw=0.8, color=col, alpha=0.28, zorder=2)
    ax.plot(T, roll_mean(y, 10), lw=1.9, color=col, label=lab, zorder=3)
style_ax(ax)
ax.set_xlabel("Time (ns)", fontsize=12.5)
ax.set_ylabel(r"Solvent accessible surface area ($\rm \AA^2$)", fontsize=12.5)
ax.set_xlim(0, T[-1])
ax.xaxis.set_major_locator(MultipleLocator(25))
ax.legend(loc="upper right", frameon=False, ncol=3)
# 界面埋藏面积 (孤立与复合物之差 = ΔSASA per component, 求和)
buried = (D["sasa_prot_iso"] - D["sasa_prot"]) + (D["sasa_dna_iso"] - D["sasa_dna"])
txt = (f"Protein  {np.nanmean(D['sasa_prot']):,.0f} ± {np.nanstd(D['sasa_prot']):,.0f} $\\rm\\AA^2$\n"
       f"ssDNA      {np.nanmean(D['sasa_dna']):,.0f} ± {np.nanstd(D['sasa_dna']):,.0f} $\\rm\\AA^2$\n"
       f"Interface (buried)  {np.nanmean(buried):,.0f} ± {np.nanstd(buried):,.0f} $\\rm\\AA^2$")
stat_box(ax, txt, loc="lower left")
axm = fig.add_subplot(gs[0, 1])
allv = np.concatenate([D["sasa_cplx"], D["sasa_prot"], D["sasa_dna"]])
lo, hi = allv.min(), allv.max()
pad = 0.06 * (hi - lo)
ys = np.linspace(lo - pad, hi + pad, 300)
for key, col in [("sasa_cplx", C_CPLX), ("sasa_prot", C_PROT), ("sasa_dna", C_DNA)]:
    kde = gaussian_kde(D[key])
    axm.plot(kde(ys), ys, color=col, lw=1.7)
    cnt, edg = np.histogram(D[key], bins=40, density=True)
    axm.barh(0.5 * (edg[:-1] + edg[1:]), cnt, height=(edg[1] - edg[0]) * 0.9,
             color=col, alpha=0.18, edgecolor="none")
axm.set_ylim(lo - pad, hi + pad)
axm.axis("off")
axm.text(0.5, 1.012, "Density", transform=axm.transAxes, ha="center",
         va="bottom", fontsize=9.5, color="#5A6670")
fig.suptitle(_ttl("溶剂可及表面积 (SASA, LCPO)"),
             fontsize=14.5, fontweight="bold", y=0.975)
save(fig, "fig04_sasa.png")


# ================================================================= FIG 5  PCA
print("[fig05] PCA")
fig = plt.figure(figsize=(12.0, 5.0))
gs = fig.add_gridspec(1, 3, width_ratios=[1.55, 1.15, 0.055], wspace=0.28,
                      left=0.062, right=0.965, top=0.855, bottom=0.145)
ax1 = fig.add_subplot(gs[0, 0])
sc = ax1.scatter(D["pc1"], D["pc2"], c=T, cmap="viridis", s=7.5,
                 alpha=0.72, edgecolors="none", zorder=3)
xmin, xmax = D["pc1"].min(), D["pc1"].max()
ymin, ymax = D["pc2"].min(), D["pc2"].max()
X, Y = np.mgrid[xmin:xmax:120j, ymin:ymax:120j]
vals = np.vstack([D["pc1"], D["pc2"]])
kde2 = gaussian_kde(vals)
Z = kde2(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
ax1.contour(X, Y, Z, levels=7, colors="#415A77", linewidths=0.85, alpha=0.65, zorder=4)
style_ax(ax1, gridy=False)
ax1.grid(ls="--", lw=0.6, color=GRIDC, alpha=0.6)
ax1.set_xlabel(f"PC1 ({var_pct[0]:.1f}%)", fontsize=12.5)
ax1.set_ylabel(f"PC2 ({var_pct[1]:.1f}%)", fontsize=12.5)
ax1.axhline(0, lw=0.9, color=GREY, alpha=0.7, zorder=1)
ax1.axvline(0, lw=0.9, color=GREY, alpha=0.7, zorder=1)
cb = fig.colorbar(sc, cax=fig.add_subplot(gs[0, 2]), pad=0.02)
cb.set_label("Time (ns)", fontsize=11)
cb.outline.set_linewidth(1.1)
cb.outline.set_edgecolor("#2B2B2B")
ax1.set_title("Conformational subspace", fontsize=12.5, pad=9)

ax2 = fig.add_subplot(gs[0, 1])
nmode = 20
idx = np.arange(1, nmode + 1)
ax2.bar(idx, var_pct[:nmode], color=C_PROT, alpha=0.85, width=0.66,
        edgecolor="white", lw=0.7, label="Individual")
cum = np.cumsum(var_pct[:nmode])
ax2.plot(idx, cum, "-o", ms=4.2, lw=1.7, color=C_ACC, label="Cumulative")
ax2.scatter(idx, cum, s=14, color=C_ACC, zorder=5)
style_ax(ax2)
ax2.set_xlabel("Principal component", fontsize=12.5)
ax2.set_ylabel("Variance explained (%)", fontsize=12.5)
ax2.set_xticks([1, 5, 10, 15, 20])
ax2.set_xlim(0.2, nmode + 0.8)
ax2.legend(loc="center right", frameon=False, fontsize=10)
ax2.text(0.97, 0.80, f"PC1–PC3 累计 {var_pct[:3].sum():.1f}%",
         transform=ax2.transAxes, fontsize=10, color="#3C4852",
         ha="right", va="center")
ax2.set_title("Eigenvalue spectrum", fontsize=12.5, pad=9)
panel_label(ax1, "A")
panel_label(ax2, "B", x=-0.03, y=1.06)
fig.suptitle(_ttl("主成分分析 (PCA, C$\\alpha$)"),
             fontsize=14.5, fontweight="bold", y=0.985)
save(fig, "fig05_pca.png")


# ================================================================= FIG 6  FEL
print("[fig06] FEL")
kT = 0.0019872 * 300.0        # kcal/mol
nb = 60
H, xe, ye = np.histogram2d(D["pc1"], D["pc2"], bins=nb, density=True)
H = gaussian_filter(H, sigma=1.6)
P = H / H.max()
with np.errstate(divide="ignore"):
    G = -kT * np.log(P)
G = np.clip(G - np.nanmin(G), 0, 6.0)
G[G >= 5.98] = np.nan
xc = 0.5 * (xe[:-1] + xe[1:])
yc = 0.5 * (ye[:-1] + ye[1:])
Xg, Yg = np.meshgrid(xc, yc)

fig = plt.figure(figsize=(12.4, 5.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.35, 0.06], wspace=0.22,
                      left=0.06, right=0.95, top=0.855, bottom=0.14)
ax1 = fig.add_subplot(gs[0, 0])
levels = np.linspace(0, 6, 13)
cf = ax1.contourf(Xg, Yg, G.T, levels=levels, cmap="Spectral_r", extend="max")
cl = ax1.contour(Xg, Yg, G.T, levels=levels[1:-1], colors="white",
                 linewidths=0.65, alpha=0.5)
style_ax(ax1, gridy=False)
ax1.set_xlabel(f"PC1 ({var_pct[0]:.1f}%)", fontsize=12.5)
ax1.set_ylabel(f"PC2 ({var_pct[1]:.1f}%)", fontsize=12.5)
# 标出全局最小与次小盆地
Gm = np.where(np.isnan(G), 1e9, G)
imin = np.unravel_index(np.argmin(Gm), Gm.shape)
ax1.plot(Xg[imin], Yg[imin], "*", ms=19, color="white", mec="#111111", mew=1.3, zorder=8)
ax1.annotate("global minimum", xy=(Xg[imin], Yg[imin]),
             xytext=(Xg[imin] + 0.09 * (Xg.max() - Xg.min()),
                     Yg[imin] + 0.10 * (Yg.max() - Yg.min())),
             fontsize=9.6, fontweight="bold", color="#111111",
             arrowprops=dict(arrowstyle="->", lw=1.2, color="#111111"))
ax1.set_title("2D free-energy landscape", fontsize=12.5, pad=9)

ax2 = fig.add_subplot(gs[0, 1], projection="3d")
Gp = np.where(np.isnan(G), 6.0, G)
surf = ax2.plot_surface(Xg, Yg, Gp.T, cmap="Spectral_r", rstride=2, cstride=2,
                        linewidth=0, antialiased=True, alpha=0.97,
                        vmin=0, vmax=6)
ax2.set_zlim(0, 6.4)
ax2.set_xlabel(f"PC1 ({var_pct[0]:.1f}%)", fontsize=10.5, labelpad=8)
ax2.set_ylabel(f"PC2 ({var_pct[1]:.1f}%)", fontsize=10.5, labelpad=8)
ax2.set_zlabel("ΔG (kcal/mol)", fontsize=10.5, labelpad=6)
ax2.view_init(elev=32, azim=-58)
ax2.xaxis.pane.fill = False
ax2.yaxis.pane.fill = False
ax2.zaxis.pane.fill = False
ax2.xaxis.pane.set_edgecolor("#DDDDDD")
ax2.yaxis.pane.set_edgecolor("#DDDDDD")
ax2.zaxis.pane.set_edgecolor("#DDDDDD")
ax2.grid(color="#E3E8EC", lw=0.6)
ax2.set_title("3D surface", fontsize=12.5, pad=4)
cb = fig.colorbar(cf, cax=fig.add_subplot(gs[0, 2]), pad=0.02)
cb.set_label("ΔG (kcal/mol)", fontsize=11.5)
cb.outline.set_linewidth(1.1)
cb.outline.set_edgecolor("#2B2B2B")
panel_label(ax1, "A")
panel_label(ax2, "B", x=-0.02, y=1.02)
fig.suptitle(_ttl("自由能形貌 (FEL)"),
             fontsize=14.5, fontweight="bold", y=0.985)
save(fig, "fig06_fel.png")


# ================================================================= FIG 7  MM-PBSA
print("[fig07] MM-PBSA")


def parse_mmpbsa(path, section):
    """返回 (GB 分量 dict, PB 分量 dict)"""
    txt = open(path, encoding="utf-8", errors="ignore").read()
    out = {}
    blocks = txt.split("Differences (Complex - Receptor - Ligand):")
    if len(blocks) < 2:
        raise RuntimeError("no delta block")
    for bi, model in [(1, "GB"), (2, "PB")]:
        if bi >= len(blocks):
            continue
        seg = blocks[bi]
        if model == "GB":
            # GB 段后面紧接 POISSON BOLTZMANN 绝对能量段, 必须切掉
            seg = seg.split("POISSON BOLTZMANN:")[0]
        d = {}
        for ln in seg.splitlines():
            p = ln.split()
            if len(p) >= 4 and re.match(r"^[A-Z0-9\-]+$", p[0]):
                try:
                    d[p[0]] = (float(p[1]), float(p[2]))
                except ValueError:
                    pass
            if "DELTA TOTAL" in ln:
                q = ln.split()
                d["DELTA TOTAL"] = (float(q[2]), float(q[3]))
        out[model] = d
    return out


def parse_decomp(path):
    """解析 per-residue 分解; 返回 dict: {model: DataFrame-like dicts}"""
    lines = open(path, encoding="utf-8", errors="ignore").read().splitlines()
    res = {}
    model = None
    i = 0
    while i < len(lines):
        ln = lines[i]
        if "Generalized Born solvent" in ln:
            model = "GB"
        elif "Poisson Boltzmann solvent" in ln:
            model = "PB"
        if ln.startswith("Total Energy Decomposition:"):
            rows = []
            j = i + 3
            while j < len(lines) and lines[j].strip():
                f = lines[j].split(",")
                if len(f) < 20:
                    j += 1
                    continue
                head = f[0].strip()
                loc = f[1].strip()
                m = re.match(r"^([A-Z0-9]+)\s+(\d+)$", head)
                if not m:
                    j += 1
                    continue
                num = np.array([float(x) for x in f[2:20]])
                rows.append(dict(name=head, loc=loc[0], num=int(m.group(2)),
                                 vdw=num[3], eel=num[6], pol=num[9],
                                 npo=num[12], tot=num[15], tot_std=num[16]))
                j += 1
            if model and rows:
                res[model] = rows
            i = j
        i += 1
    return res


_MB_PATH = os.path.join(RES, "FINAL_RESULTS.dat")
_DC_PATH = os.path.join(RES, "FINAL_DECOMP_MMPBSA.dat")
if os.path.isfile(_MB_PATH) and os.path.isfile(_DC_PATH):
    MB = parse_mmpbsa(_MB_PATH, "d")
    DC = parse_decomp(_DC_PATH)
    gb_rows = DC.get("GB", [])
    print(f"   MM-PBSA GB ΔG={MB['GB']['DELTA TOTAL'][0]:.2f}  "
          f"PB ΔG={MB['PB']['DELTA TOTAL'][0]:.2f}  decomp rows={len(gb_rows)}")
else:
    MB = DC = None
    gb_rows = []
    print("[跳过 fig07] 缺 MM-PBSA 产物（FINAL_RESULTS.dat / FINAL_DECOMP_MMPBSA.dat）")

if MB is None:
    # 缺 MMPBSA 数据: 输出占位说明图, 保 fig 编号连续（fig21/24 同理由 plot2/3 处理）
    fig, ax = plt.subplots(figsize=(8.4, 3.6))
    ax.axis("off")
    ax.text(0.5, 0.5, "fig07 skipped — no MM-PBSA data\n(run stage5 or set DO_MMPBSA=yes)",
            ha="center", va="center", fontsize=13, color="#8A8F98")
    save(fig, "fig07_mmpbsa.png")
    print("[fig07] placeholder saved")
else:
    # ---- fig07 原绘图体（仅 MB 可用时执行）----

    fig = plt.figure(figsize=(12.6, 5.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.0], wspace=0.235,
                          left=0.058, right=0.975, top=0.855, bottom=0.135)

    # --- Panel A: 能量项 (GB vs PB 分组柱), PB 非极性 = ENPOLAR + EDISPER
    ax1 = fig.add_subplot(gs[0, 0])
    pb_np = tuple(np.sum([MB["PB"].get(k, (0.0, 0.0))[i] for k in ("ENPOLAR", "EDISPER")])
                  for i in (0, 1))
    terms = [("VDWAALS", "ΔE$_{vdW}$"), ("EEL", "ΔE$_{elec}$"),
             ("EGB", "ΔG$_{polar}^{GB}$"), ("EPB", "ΔG$_{polar}^{PB}$"),
             ("ESURF", "ΔG$_{np}^{GB}$"), (None, "ΔG$_{np}^{PB}$"),
             ("DELTA TOTAL", "ΔG$_{bind}$")]
    gbv, gbe, pbv, pbe, labels, cols = [], [], [], [], [], []
    cmap = [C_DNA, C_DNA, C_PROT, C_PROT, C_CPLX, C_CPLX, C_ACC]
    for (k, lab), c in zip(terms, cmap):
        a = MB["GB"].get(k) if k else None
        b = pb_np if k is None else MB["PB"].get(k)
        if a is None and b is None:
            continue
        labels.append(lab)
        gbv.append(a[0] if a else np.nan)
        gbe.append(a[1] if a else 0.0)
        pbv.append(b[0] if b else np.nan)
        pbe.append(b[1] if b else 0.0)
        cols.append(c)
    x = np.arange(len(labels))
    w = 0.37
    ax1.bar(x - w / 2, gbv, w, yerr=gbe, color=C_PROT, alpha=0.92, label="GB (igb=5)",
            edgecolor="white", lw=0.9, error_kw=dict(ecolor="#2B2B2B", capsize=3.5, lw=1.1))
    ax1.bar(x + w / 2, pbv, w, yerr=pbe, color=C_ACC, alpha=0.92, label="PB (istrng=0.15)",
            edgecolor="white", lw=0.9, error_kw=dict(ecolor="#2B2B2B", capsize=3.5, lw=1.1))
    style_ax(ax1)
    ax1.axhline(0, lw=1.15, color=INK)
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, fontsize=10.2)
    ax1.set_ylabel("Energy (kcal/mol)", fontsize=12.5)
    ax1.legend(loc="lower left", frameon=False, fontsize=10.5)
    # ΔG 总量标注
    ax1.text(0.985, 0.965,
             f"ΔG$_{{bind}}$(GB) = {MB['GB']['DELTA TOTAL'][0]:.1f} ± {MB['GB']['DELTA TOTAL'][1]:.1f}\n"
             f"ΔG$_{{bind}}$(PB) = {MB['PB']['DELTA TOTAL'][0]:.1f} ± {MB['PB']['DELTA TOTAL'][1]:.1f}  kcal/mol",
             transform=ax1.transAxes, ha="right", va="top", fontsize=10.6,
             fontweight="bold", color="#1A1A1A",
             bbox=dict(boxstyle="round,pad=0.42", fc="#FFF7E6", ec=C_ACC, lw=1.2))
    ax1.set_title("Binding free-energy components", fontsize=12.5, pad=9)

    # --- Panel B: per-residue 分解 (Top 稳定化残基)
    ax2 = fig.add_subplot(gs[0, 1])
    rows = [r for r in gb_rows]
    rows.sort(key=lambda r: r["tot"])
    top = rows[:22][::-1]
    ylab, ytot, ystd, ycol = [], [], [], []
    for r in top:
        is_dna = r["loc"] == "L"
        if is_dna:
            pos = r["num"] - 254
            ylab.append(f"dN{pos}")
            ycol.append(C_DNA)
        else:
            ylab.append(f"{r['name'].split()[0].title()}{r['num'] + OFF}")
            ycol.append(C_PROT if r["tot"] < -1.0 else "#8FB6DC")
        ytot.append(r["tot"])
        ystd.append(r["tot_std"])
    y = np.arange(len(top))
    ax2.barh(y, ytot, xerr=ystd, color=ycol, height=0.72, edgecolor="white", lw=0.7,
             error_kw=dict(ecolor="#2B2B2B", capsize=2.6, lw=0.9))
    style_ax(ax2, gridy=False, gridx=True)
    ax2.axvline(0, lw=1.15, color=INK)
    ax2.set_yticks(y)
    ax2.set_yticklabels(ylab, fontsize=9.2)
    ax2.set_xlabel("ΔG$_{residue}$ (kcal/mol, GB idecomp=1)", fontsize=12)
    ax2.set_title("Per-residue decomposition (Top 22 stabilizing)", fontsize=12.5, pad=9)
    ax2.text(0.98, 0.025, "blue = protein   ·   orange = ssDNA nucleotide",
             transform=ax2.transAxes, fontsize=9, color="#5A6670", ha="right")
    panel_label(ax1, "A")
    panel_label(ax2, "B", x=-0.025, y=1.06)
    fig.suptitle(_ttl("MM-PBSA 结合自由能 (GB igb=5 / PB istrng=0.15)"),
                 fontsize=14.5, fontweight="bold", y=0.985)
    save(fig, "fig07_mmpbsa.png")


# ================================================================= FIG 8  Mn 配位
print("[fig08] Mn coordination")
_mn_keys = ("mn_h71", "mn_h75", "mn_e102a", "mn_e102b")
if not all(k in D for k in _mn_keys):
    fig, ax = plt.subplots(figsize=(8.4, 3.6))
    ax.axis("off")
    ax.text(0.5, 0.5, "fig08 skipped — no Mn$^{2+}$ coordination data\n(metal-free system, or no Mn in input PDB)",
            ha="center", va="center", fontsize=13, color="#8A8F98")
    save(fig, "fig08_mn_coordination.png")
    print("[fig08] placeholder saved")
else:
    fig = plt.figure(figsize=(10.2, 5.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[5.4, 1.15], wspace=0.07,
                          left=0.072, right=0.972, top=0.855, bottom=0.145)
    ax = fig.add_subplot(gs[0, 0])
    series = [("mn_h71", "His92 N$_{\\epsilon2}$", C_PROT),
              ("mn_h75", "His96 N$_{\\epsilon2}$", C_PUR),
              ("mn_e102a", "Glu123 O$_{\\epsilon1}$", C_DNA),
              ("mn_e102b", "Glu123 O$_{\\epsilon2}$", C_ACC)]
    for key, lab, col in series:
        y = D[key]
        ax.plot(T, y, lw=0.8, color=col, alpha=0.28, zorder=2)
        ax.plot(T, roll_mean(y, 10), lw=1.85, color=col, label=lab, zorder=3)
    ax.axhspan(1.8, 2.4, color="#B7E4C7", alpha=0.42, zorder=0)
    ax.axhline(2.10, ls="--", lw=1.1, color="#2D6A4F", alpha=0.8, zorder=1)
    style_ax(ax)
    ax.set_xlabel("Time (ns)", fontsize=12.5)
    ax.set_ylabel(r"Mn$^{2+}$—ligand distance ($\rm \AA$)", fontsize=12.5)
    ax.set_xlim(0, T[-1])
    ax.set_ylim(1.4, 6.2)
    ax.xaxis.set_major_locator(MultipleLocator(25))
    ax.yaxis.set_major_locator(MultipleLocator(1.0))
    ax.legend(loc="upper right", frameon=False, ncol=2, fontsize=9.8)
    lines = []
    for key, lab, col in series:
        y = D[key]
        lines.append(f"{lab:<18s} {np.nanmean(y):.2f} ± {np.nanstd(y):.2f} $\\rm\\AA$")
    stat_box(ax, "\n".join(lines), loc="upper left", fs=8.8)

    axm = fig.add_subplot(gs[0, 1])
    allv = np.concatenate([D[k] for k, _, _ in series])
    lo, hi = 1.6, min(5.0, np.nanpercentile(allv, 99.5))
    ys = np.linspace(lo, hi, 300)
    for key, lab, col in series:
        d = D[key][(D[key] >= lo) & (D[key] <= hi)]
        if len(d) < 10:
            continue
        kde = gaussian_kde(d)
        axm.plot(kde(ys), ys, color=col, lw=1.7)
        cnt, edg = np.histogram(d, bins=45, density=True)
        axm.barh(0.5 * (edg[:-1] + edg[1:]), cnt, height=(edg[1] - edg[0]) * 0.9,
                 color=col, alpha=0.18, edgecolor="none")
    axm.set_ylim(lo, hi)
    axm.axis("off")
    axm.text(0.5, 1.012, "Density", transform=axm.transAxes, ha="center",
             va="bottom", fontsize=9.5, color="#5A6670")
    fig.suptitle(_ttl("Mn$^{2+}$ 三齿配位稳定性"),
                 fontsize=14.5, fontweight="bold", y=0.975)
    save(fig, "fig08_mn_coordination.png")

print("\n[done] all figures ->", OUT)
