# -*- coding: utf-8 -*-
"""plot_pretty_figs3.py — 论文级 MD 后分析图集 III (fig23-30)

复用 plot_pretty_figs.py (fig01-08) 样式系统。运行:
    python plot_pretty_figs3.py
图集 III 覆盖"论文高频"主题:
  fig23 DCCM 蛋白运动互相关矩阵 (CHPC matrix correl, rms-fit 后)
  fig24 MM-PBSA 蛋白 per-residue ΔG 分解 (FINAL_DECOMP 的 R 侧)
  fig25 界面 H 键 occupancy 分布 (直方图 + 累积)
  fig26 界面接触 (native/non-native) + 总 H 键 双通道时序
  fig27 DSSP 二级结构时间演化栅栏 (1500 帧 x 254 残基)
  fig28 RMSF vs 晶体 B-factor 校验 (8SLN CA, 当量 RMSF)
  fig29 ssDNA 形态三指标: 端距 / Rg / 弯曲角
  fig30 Mn2+ 配位 RDF + 相对配位累积
"""
import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import plot_pretty_figs as P

C_PROT, C_DNA, C_CPLX = P.C_PROT, P.C_DNA, P.C_CPLX
C_ACC, C_PUR, C_TEAL = P.C_ACC, P.C_PUR, P.C_TEAL
GREY, INK = P.GREY, P.INK
RES, RAW, OUT = P.RES, os.path.join(P.RES, "cpptraj_raw"), P.OUT
save = P.save
style_ax = P.style_ax
panel_label = P.panel_label
stat_box = P.stat_box
roll_mean = P.roll_mean
KEY = P.KEY_RES
HELIX = P.HELIX_CAT        # 催化螺旋 (62,80)
LOOP = P.LOOP_MOD          # 无序环 (168,181)
S1 = "TCATGAGCAGTTTTTTGTTTTTTT"
NS = 1500
DT = 0.1                   # ns/帧
T = np.arange(NS) * DT


def read(path, skip=1):
    return np.loadtxt(path, skiprows=skip)


def mstd(y, nd=2):
    return f"{np.mean(y):.{nd}f} ± {np.std(y):.{nd}f}"


# ================================================================ FIG 23 DCCM
print("[fig23] DCCM protein cross-correlation")
DCCM = read(os.path.join(RAW, "dccm_prot.dat"), skip=0)
n = DCCM.shape[0]
assert n == 254, f"DCCM shape {DCCM.shape}"
# 对角线 1; 截断极端值做对称色标
off = DCCM[~np.eye(n, dtype=bool)]
vlim = max(abs(np.percentile(off, 1)), abs(np.percentile(off, 99))) * 0.99
# 功能区平均相关 (跨域耦合信号)
blk_cat = DCCM[HELIX[0]-1:HELIX[1], 180:254].mean()   # 催化螺旋 x C端核酸结合
blk_n = DCCM[HELIX[0]-1:HELIX[1], 55:150].mean()       # 催化螺旋 x 核心1区
fig = plt.figure(figsize=(8.6, 7.4))
gs = fig.add_gridspec(1, 1, left=0.11, right=0.93, top=0.90, bottom=0.09)
ax = fig.add_subplot(gs[0])
im = ax.imshow(DCCM, cmap="coolwarm", vmin=-vlim, vmax=vlim,
               interpolation="nearest", aspect="auto")
# 功能单元方块: 催化螺旋 / 核心1区 / C端结合区
for (a, b, lab, col) in [(HELIX[0], HELIX[1], "催化螺旋 HEXXH", C_ACC),
                         (55, 155, "核心1 (GAF)", C_PROT),
                         (180, 254, "核心2/Patch1 (C)", C_PUR)]:
    ax.add_patch(plt.Rectangle((a-0.5, a-0.5), b-a+1, b-a+1, fill=False,
                               edgecolor=col, lw=1.6, ls="--", alpha=0.85))
    ax.text(b + 2.5, a - 0.5, lab, color=col, fontsize=9, va="center",
            fontweight="bold", rotation=90)
# 关键残基刻度
ticks = sorted(set(list(range(1, 255, 50)) + [64, 67, 106, 146, 186, 199, 232, 246, 254]))
labels = [f"{t}\n{KEY[t][0]}" if t in KEY else str(t) for t in ticks]
ax.set_xticks(ticks); ax.set_yticks(ticks)
ax.set_xticklabels(labels, fontsize=8)
ax.set_yticklabels(labels, fontsize=8)
ax.set_xlabel("蛋白残基 (seq)", fontsize=12)
ax.set_ylabel("蛋白残基 (seq)", fontsize=12)
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
cb.set_label("运动互相关 C(i,j)", fontsize=11)
ax.set_title("DCCM：蛋白 CA 运动互相关矩阵（rms-fit 去除刚体后）", fontsize=13.5, pad=8)
stat_box(ax, f"⟨C⟩$_{{催化×C端}}$ = {blk_cat:+.2f}\n"
             f"⟨C⟩$_{{催化×核心1}}$ = {blk_n:+.2f}",
         loc="lower right")
fig.suptitle("PprI(WT)·S1-ssDNA·Mn$^{2+}$   150 ns MD", fontsize=14, y=0.985)
save(fig, "fig23_dccm.png")
print(f"   DCCM vlim ±{vlim:.2f}; cat×Cterm {blk_cat:+.2f}; cat×core1 {blk_n:+.2f}")

# ================================================================ FIG 24 蛋白 per-residue 分解
print("[fig24] MM-PBSA protein per-residue")
gb_rows = P.parse_decomp(os.path.join(RES, "FINAL_DECOMP_MMPBSA.dat")).get("GB", [])
prot = [r for r in gb_rows if r.get("loc") == "R"]
prot.sort(key=lambda r: r["num"])
seqs = np.array([r["num"] for r in prot])
tots = np.array([r["tot"] for r in prot])
key_ids = [k for k in KEY if k in seqs]
fig = plt.figure(figsize=(12.4, 5.2))
ax = fig.add_subplot(111)
ax.axvspan(0.5, 155.5, color=C_PROT, alpha=0.045, lw=0)
ax.axvspan(155.5, 254.5, color=C_PUR, alpha=0.045, lw=0)
ax.text(78, ax.get_ylim()[1] if False else 1.02, "", fontsize=1)
cols = np.where(tots < 0, C_PROT, "#E8B4B8")
ax.bar(seqs, tots, color=cols, width=0.85, alpha=0.88)
ax.axhline(0, lw=1.2, color=INK)
for s in seqs[::25]:
    pass
# 关键残基竖线
for k in key_ids:
    ax.axvline(k, color=GREY, ls=":", lw=0.7, alpha=0.8)
    ax.text(k, ax.get_ylim()[0] - (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.03,
            KEY[k][0], rotation=60, fontsize=7.2, color="#5A6670",
            ha="center", va="top")
ax.text(78, 1.0, "GAF 域", transform=ax.get_xaxis_transform(),
        ha="center", color=C_PROT, fontsize=10, fontweight="bold", alpha=0.75)
ax.text(205, 1.0, "C 端蛋白酶域", transform=ax.get_xaxis_transform(),
        ha="center", color=C_PUR, fontsize=10, fontweight="bold", alpha=0.75)
ax.set_xlim(0.5, 254.5)
ax.set_xlabel("蛋白残基 (seq)", fontsize=12.5)
ax.set_ylabel("ΔG$_{res}$ (kcal/mol, GB)", fontsize=12.5)
style_ax(ax)
top = sorted(zip(seqs, tots), key=lambda z: z[1])[:6]
lab = "  ".join(f"{KEY.get(s,(str(s),''))[0]} {v:.1f}" for s, v in top)
stat_box(ax, f"ΣΔG = {tots.sum():.1f} kcal/mol\n稳定残基: {lab}", loc="upper left")
ax.set_title("MM-PBSA 蛋白 per-residue 分解（负值 = 稳定结合贡献）", fontsize=13)
save(fig, "fig24_mmpbsa_prot_residue.png")
print(f"   ΣΔG {tots.sum():.1f}; top {[(KEY.get(int(s),('',))[0] or int(s), round(float(v),1)) for s,v in top]}")

# ================================================================ FIG 25 H 键 occupancy 分布
print("[fig25] H-bond occupancy distribution")
rows = []
for ln in open(os.path.join(RAW, "hb_all_avg.dat"), encoding="utf-8", errors="ignore"):
    if ln.startswith("#") or not ln.strip():
        continue
    f = ln.split()
    if len(f) < 7:
        continue
    acc, don = f[0], f[2]
    m_acc = re.match(r"[A-Z0-9]+_(\d+)@.+", acc)
    m_don = re.match(r"[A-Z0-9]+_(\d+)@.+", don)
    if not m_acc or not m_don:
        continue
    ra, rd = int(m_acc.group(1)), int(m_don.group(1))
    if (ra <= 254 < rd) or (rd <= 254 < ra):       # 界面
        rows.append((float(f[4]), float(f[5])))     # frac, dist
rows.sort(reverse=True)
occ = np.array([r[0] * 100 for r in rows])
dist = np.array([r[1] for r in rows])
fig = plt.figure(figsize=(9.6, 5.2))
gs = fig.add_gridspec(1, 2, width_ratios=[1.7, 1.0], wspace=0.13,
                      left=0.10, right=0.98, top=0.86, bottom=0.14)
ax = fig.add_subplot(gs[0])
bins = np.arange(0, 105, 5)
ax.hist(occ, bins=bins, color=C_PROT, alpha=0.75, edgecolor="white", lw=0.6,
        zorder=3)
cum = np.cumsum(np.histogram(occ, bins=bins)[0]) / len(occ) * 100
ax2 = ax.twinx()
ax2.plot(bins[:-1] + 2.5, cum, color=C_ACC, lw=2.2, marker="o", ms=3, zorder=4)
ax2.set_ylim(0, 105)
ax2.set_ylabel("累积占比 (%)", fontsize=11, color=C_ACC)
ax2.spines["top"].set_visible(False)
ax.set_xlabel("界面 H 键 occupancy (%, 300 帧)", fontsize=12)
ax.set_ylabel("氢键数目", fontsize=12)
style_ax(ax, gridy=True)
for th, col, nm in [(20, "#C9A2C8", "occ<20%"), (60, "#7B4EA8", "occ≥60%")]:
    n = (occ >= th if th > 20 else occ < th).sum()
    ax.axvline(th, ls="--", lw=0.9, color=col)
    ax.text(th, ax.get_ylim()[1] * 0.96, f" {nm}: n={n}", fontsize=8.5,
            color=col, fontweight="bold")
# 右侧: occupancy vs distance 2D 关系 (散射密度)
axr = fig.add_subplot(gs[1])
hb = axr.hexbin(dist, occ, gridsize=34, cmap="Blues", mincnt=1, extent=(2.4, 4.2, 0, 100))
axr.axhline(50, ls="--", lw=0.8, color=GREY)
cb = fig.colorbar(hb, ax=axr, fraction=0.05, pad=0.03)
cb.set_label("n")
axr.set_xlabel("Avg heavy-atom dist (Å)", fontsize=11)
axr.set_ylabel("Occupancy (%)", fontsize=11)
axr.set_title("occ × dist", fontsize=11)
strong = ((dist < 3.2) & (occ >= 50)).sum()
ax.set_title(f"界面 H 键 occupancy 分布（n={len(occ)}；强且持久 <3.2Å & ≥50%: {strong}）",
             fontsize=13)
save(fig, "fig25_hbond_occupancy_dist.png")
print(f"   n={len(occ)}; occ≥60%: {(occ>=60).sum()}; <20%: {(occ<20).sum()}")

# ================================================================ FIG 26 双通道时序
print("[fig26] contacts + Hbond dual time series")
nat = read(os.path.join(RAW, "nat_series.dat"))          # 300 帧
hb = read(os.path.join(RAW, "hb_all_series.dat"))         # 300 帧
tnat = (nat[:, 0] - 1) * 5 * DT                          # stride 5
thb = (hb[:, 0] - 1) * 5 * DT
fig = plt.figure(figsize=(11.4, 5.6))
gs = fig.add_gridspec(1, 1, left=0.09, right=0.985, top=0.87, bottom=0.13)
ax = fig.add_subplot(gs[0])
ax.plot(tnat, nat[:, 1], color=C_PROT, lw=1.0, alpha=0.35, zorder=2)
ax.plot(tnat, roll_mean(nat[:, 1], 15), color=C_PROT, lw=2.4, zorder=3,
        label="native contacts")
ax.plot(tnat, nat[:, 2], color="#B7B7B7", lw=0.8, alpha=0.5, zorder=1)
ax.plot(tnat, roll_mean(nat[:, 2], 15), color="#7A8B99", lw=1.8, ls="--", zorder=2,
        label="non-native contacts")
ax.set_xlabel("Time (ns)", fontsize=12.5)
ax.set_ylabel("原子接触数 (nativecontacts)", fontsize=12, color=C_PROT)
ax.tick_params(axis="y", labelcolor=C_PROT)
ax2 = ax.twinx()
ax2.plot(thb, hb[:, 1], color=C_DNA, lw=1.0, alpha=0.3, zorder=2)
ax2.plot(thb, roll_mean(hb[:, 1], 15), color=C_DNA, lw=2.2, zorder=3,
         label="H-bonds (all)")
ax2.set_ylabel("界面 H 键总数", fontsize=12, color=C_DNA)
ax2.tick_params(axis="y", labelcolor=C_DNA)
style_ax(ax)
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=10, frameon=False)
ax.set_title("界面接触与氢键双通道动态（滚动均值 7.5 ns）", fontsize=13.5)
stat_box(ax, f"⟨native⟩ = {nat[:,1].mean():.0f} ± {nat[:,1].std():.0f}\n"
             f"⟨H-bond⟩ = {hb[:,1].mean():.0f} ± {hb[:,1].std():.0f}",
         loc="upper left")
save(fig, "fig26_contacts_hbond_dual.png")
print(f"   native {nat[:,1].mean():.0f}±{nat[:,1].std():.0f}; hb {hb[:,1].mean():.0f}±{hb[:,1].std():.0f}")

# ================================================================ FIG 27 DSSP timeline
print("[fig27] DSSP secondary-structure timeline")
d = np.loadtxt(os.path.join(RAW, "dssp.dat"), skiprows=1)
if d.ndim == 2 and d.shape[1] == 255:
    d = d[:, 1:]                       # strip Frame col
SS = d.astype(int)                     # 1500 x 254
code2ss = {0: "coil", 1: "Extended", 2: "Bridge", 3: "3-10",
           4: "Alpha", 5: "Pi", 6: "Turn", 7: "Bend"}
cmap = ListedColormap(["#F1F3F5", "#E8C23A", "#C9971C", "#8FB8E8",
                       "#D64541", "#B071D8", "#4BAE5E", "#5B9BD5"])
norm = BoundaryNorm(np.arange(-0.5, 8.5, 1), cmap.N)
step = 10
S = SS[::step]                         # 150 列
fig = plt.figure(figsize=(12.6, 6.6))
ax = fig.add_subplot(111)
im = ax.imshow(S.T, aspect="auto", cmap=cmap, norm=norm,
               extent=[0, S.shape[0] * step * DT, 0.5, 254.5],
               interpolation="nearest")
ax.axhspan(HELIX[0] - 0.5, HELIX[1] + 0.5, color=C_ACC, alpha=0.10, lw=0)
ax.axhspan(LOOP[0] - 0.5, LOOP[1] + 0.5, color=GREY, alpha=0.16, lw=0)
ax.text(S.shape[0] * step * DT * 0.995, (HELIX[0] + HELIX[1]) / 2,
        "催化螺旋 HEXXH", fontsize=9, color="#B7791F", fontweight="bold",
        ha="right", va="center", rotation=90)
ax.text(S.shape[0] * step * DT * 0.995, (LOOP[0] + LOOP[1]) / 2,
        "无序环 (建模)", fontsize=8.5, color="#5A6670",
        ha="right", va="center", rotation=90)
for k, (nm, grp) in KEY.items():
    col = P.GROUP_COLOR.get(grp, GREY)
    ax.axhline(k, color=col, ls=":", lw=0.5, alpha=0.55)
ax.set_yticks(list(range(1, 255, 20)) + sorted(k for k in KEY if k % 20))
ax.set_yticklabels([str(t) for t in ax.get_yticks()], fontsize=8)
ax.set_ylim(254.5, 0.5)
ax.set_xlabel("Time (ns)", fontsize=12.5)
ax.set_ylabel("蛋白残基 (seq)", fontsize=12.5)
cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, ticks=[0, 1, 2, 3, 4, 5, 6, 7])
cbar.ax.set_yticklabels([code2ss[i] for i in range(8)], fontsize=8)
# 每帧 alpha 占比
alpha_frac = (SS == 4).mean(axis=1)
ax2 = ax.twinx()
tt = np.arange(len(alpha_frac)) * DT
ax2.plot(tt[::step], alpha_frac[::step], color="#2C3E50", lw=1.2, alpha=0.8)
ax2.set_ylim(0, 0.6)
ax2.set_ylabel("α-helix 占比", fontsize=10.5)
ax2.tick_params(axis="y", labelcolor="#2C3E50", labelsize=8)
ax2.spines["top"].set_visible(False)
ax.set_title("DSSP 二级结构时间演化栅栏（1500 帧；彩色 = 蛋白关键位点）", fontsize=13.5)
save(fig, "fig27_dssp_timeline.png")
print(f"   SS {SS.shape}; alpha_frac {alpha_frac.mean():.2f}±{alpha_frac.std():.2f}")

# ================================================================ FIG 28 RMSF vs B-factor
print("[fig28] RMSF vs crystal B-factor")
rmsf = read(os.path.join(RES, "rmsf_prot.dat"))
bf = np.loadtxt(os.path.join(RES, "bfactor_prot.dat"), skiprows=1)   # seq pdb B (nan)
B = np.full(254, np.nan)
B[bf[:, 0].astype(int) - 1] = bf[:, 2]           # seq -> index
B_eq = np.sqrt(3.0 * B / (8.0 * np.pi ** 2))          # 当量 RMSF (Å)
x = np.arange(1, 255)
good = ~np.isnan(B_eq)
fig = plt.figure(figsize=(12.4, 5.0))
ax = fig.add_subplot(111)
ax.fill_between(x, rmsf[:, 1], color=C_PROT, alpha=0.28, lw=0)
ax.plot(x, rmsf[:, 1], color=C_PROT, lw=1.6, label="MD RMSF (CA)")
ax.plot(x[good], B_eq[good], color=C_DNA, lw=1.7, alpha=0.9,
        label="crystal B-factor 当量 RMSF")
r = np.corrcoef(rmsf[good, 1], B_eq[good])[0, 1]
ax.axvspan(LOOP[0], LOOP[1], color=GREY, alpha=0.18, lw=0)
ax.text((LOOP[0] + LOOP[1]) / 2, ax.get_ylim()[1] * 0.94, "无序环\n(建模)",
        ha="center", fontsize=8.5, color="#5A6670")
for k, (nm, grp) in KEY.items():
    col = P.GROUP_COLOR.get(grp, GREY)
    ax.axvline(k, color=col, ls=":", lw=0.8, alpha=0.75)
    ax.text(k, ax.get_ylim()[0] + 0.04, KEY[k][0], rotation=55, fontsize=7.2,
            color=col, ha="right", va="bottom", fontweight="bold")
ax.set_xlim(1, 254)
ax.set_xlabel("蛋白残基 (seq)", fontsize=12.5)
ax.set_ylabel("RMSF / 当量 RMSF (Å)", fontsize=12.5)
style_ax(ax)
ax.legend(loc="upper left", fontsize=10, frameon=False)
stat_box(ax, f"Pearson r = {r:.2f}  (n={good.sum()})", loc="upper right")
ax.set_title("MD 柔性 vs 晶体 B-factor 校验（B → √(3B/8π²) 当量化）", fontsize=13)
save(fig, "fig28_rmsf_vs_bfactor.png")
print(f"   Pearson r = {r:.2f} (n={good.sum()})")

# ================================================================ FIG 29 DNA 形态三指标
print("[fig29] ssDNA morphology: EE / Rg / bend")
ee = read(os.path.join(RAW, "dna_ee.dat"))[:, 1]
rg = read(os.path.join(RES, "rg_dna.dat"))[:, 1]
bd = read(os.path.join(RAW, "dna_bend.dat"))[:, 1]
fig = plt.figure(figsize=(11.8, 7.6))
gs = fig.add_gridspec(3, 1, hspace=0.36, left=0.09, right=0.985,
                      top=0.94, bottom=0.07)
pan = [(ee, C_TEAL, "端到端距离", "Å"),
       (rg, C_DNA, "回旋半径 R$_g$", "Å"),
       (bd, C_PUR, "弯曲角 (5′-中-3′)", "°")]
for i, (y, col, nm, unit) in enumerate(pan):
    ax = fig.add_subplot(gs[i])
    ax.plot(T, y, lw=0.9, color=col, alpha=0.4)
    ax.plot(T, roll_mean(y, 30), lw=2.2, color=col)
    if i == 2:
        ax.axhline(180, ls="--", lw=0.8, color=GREY)
        ax.text(1, 181.5, "直线 (180°)", fontsize=8.5, color="#5A6670")
    ax.set_ylabel(f"{nm} ({unit})", fontsize=11.5, color=col)
    panel_label(ax, chr(97 + i))
    stat_box(ax, mstd(y), loc="upper right", fs=9.5)
    if i < 2:
        ax.set_xticklabels([])
    else:
        ax.set_xlabel("Time (ns)", fontsize=12.5)
    style_ax(ax)
fig.suptitle("ssDNA 形态动态：端距 / 回旋半径 / 弯曲角（滚动均值 3 ns）",
             fontsize=14, y=0.985)
save(fig, "fig29_dna_morphology.png")
print(f"   EE {mstd(ee)}; Rg {mstd(rg)}; bend {mstd(bd,0)}")

# ================================================================ FIG 30 Mn RDF
print("[fig30] Mn2+ RDF")
rdf = read(os.path.join(RES, "mn_rdf.dat"), skip=1)
r_ = rdf[:, 0]
g = rdf[:, 1]
fig = plt.figure(figsize=(8.2, 5.2))
ax = fig.add_subplot(111)
ax.plot(r_, g, color=C_ACC, lw=2.2)
ax.fill_between(r_, g, color=C_ACC, alpha=0.22)
# 第一配位壳层峰
peak_i = np.argmax(g[(r_ >= 1.5) & (r_ <= 3.5)])
r_peak = r_[(r_ >= 1.5) & (r_ <= 3.5)][peak_i]
ax.axvline(r_peak, ls="--", lw=1.1, color=C_ACC)
ax.text(r_peak + 0.15, g.max() * 0.9, f"第一壳层 {r_peak:.1f} Å",
        fontsize=10.5, color="#B7791F", fontweight="bold")
# 晶体配位键参考 (fig08: H71/H75/E102 ~2.2/2.25/2.77)
for xr, lab in [(2.23, "H71"), (2.25, "H75"), (2.77, "E102")]:
    ax.axvline(xr, ls=":", lw=0.9, color=GREY, alpha=0.8)
    ax.text(xr, g.max() * 0.06, lab, fontsize=8, color="#5A6670",
            ha="center", rotation=90)
ax.set_xlim(0, 10)
ax.set_xlabel("r (Å)", fontsize=12.5)
ax.set_ylabel("g(r)  Mn$^{2+}$ → 蛋白原子", fontsize=12.5)
style_ax(ax)
ax.set_title("Mn$^{2+}$ 配位环境径向分布（RDF, 150 ns）", fontsize=13)
save(fig, "fig30_mn_rdf.png")
print(f"   first-shell peak {r_peak:.2f} Å; gmax {g.max():.2f}")

print("ALL DONE (fig23-30)")
