#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PprI(WT)·S1_ssDNA·Mn 150 ns MD —— 图集扩展 fig09..fig22 (14 张)

import plot_pretty_figs 以复用其全部样式/配色/残基注释与数据 (D 字典)。
新数据源: cpptraj_raw/ (dssp, phi/psi, nativecontacts, hbond, pocket, dna_ee) + prod.out
"""
import os, re
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.ticker import MultipleLocator, FuncFormatter

import plot_pretty_figs as P          # 重跑 fig01-08 (同款样式) 并复用其命名空间
from plot_pretty_figs import (C_PROT, C_DNA, C_CPLX, C_ACC, C_PUR, C_TEAL,
                              GREY, INK, GRIDC, PALETTE, KEY_RES, OFF,
                              HELIX_CAT, LOOP_MOD, GROUP_COLOR,
                              read_dat, roll_mean, style_ax, stat_box,
                              save, RES, OUT, D, NFR, DT_NS, T,
                              parse_decomp)

RAW = os.path.join(RES, "cpptraj_raw")
S1 = "TCATGAGCAGTTTTTTGTTTTTTT"       # 24 nt 靶序列
print(f"[fig09-22] loaded base figs, {NFR} frames, {T[-1]:.0f} ns")

# ================================================================ FIG 09 能量收敛
print("[fig09] energy convergence")
def parse_prodout(path):
    txt = open(path, encoding="utf-8", errors="ignore").read().splitlines()
    nstep_pos = []
    for i, ln in enumerate(txt):
        if re.search(r"NSTEP\s*=\s*\d+\s+TIME\(PS\)", ln):
            nstep_pos.append(i)
    recs = []
    for idx, i in enumerate(nstep_pos):
        m = re.search(r"NSTEP\s*=\s*(\d+)\s+TIME\(PS\)\s*=\s*([\d.]+)\s+TEMP\(K\)\s*=\s*([\d.]+)\s+PRESS\s*=\s*([-\d.]+)", txt[i])
        r = dict(nstep=int(m.group(1)), t=float(m.group(2)),
                 temp=float(m.group(3)), press=float(m.group(4)))
        jend = nstep_pos[idx + 1] if idx + 1 < len(nstep_pos) else min(i + 16, len(txt))
        for k in range(i + 1, min(jend, i + 16)):
            e = re.search(r"Etot\s*=\s*([-\d.]+)\s+EKtot\s*=\s*([-\d.]+)\s+EPtot\s*=\s*([-\d.]+)", txt[k])
            if e and "etot" not in r:
                r.update(etot=float(e.group(1)), ek=float(e.group(2)), ep=float(e.group(3)))
            v = re.search(r"VOLUME\s*=\s*([\d.]+)", txt[k])
            if v and "vol" not in r:
                r["vol"] = float(v.group(1))
            d = re.search(r"Density\s*=\s*([\d.]+)", txt[k])
            if d and "dens" not in r:
                r["dens"] = float(d.group(1))
        if "etot" in r and "vol" in r and "dens" in r:
            recs.append(r)
    return recs

ER = parse_prodout(os.path.join(RES, "prod.out"))
print(f"   {len(ER)} energy blocks; t {ER[0]['t']:.0f}-{ER[-1]['t']:.0f} ps")
te = np.array([r["t"] for r in ER]) / 1000.0      # ns
temp, press, dens = ([np.array([r[k] for r in ER]) for k in ("temp", "press", "dens")])
etot, ek, ep = ([np.array([r[k] for r in ER]) for k in ("etot", "ek", "ep")])
fig = plt.figure(figsize=(11.6, 6.6))
gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.30, left=0.075, right=0.985, top=0.875, bottom=0.10)
titles = [("TEMP(K)", temp, C_DNA, "Temperature"), ("PRESS(bar)", press, C_PUR, "Pressure"),
          ("Density (g/cm³)", dens, C_CPLX, "Density"),
          ("E_total (kcal/mol)", etot, C_PROT, "Total energy"),
          ("E_potential", ep, C_TEAL, "Potential energy"), ("E_kinetic", ek, C_ACC, "Kinetic energy")]
for a, (lab, yy, col, cap) in enumerate(titles):
    ax = fig.add_subplot(gs[a // 3, a % 3])
    ax.plot(te, yy, lw=0.9, color=col, alpha=0.55)
    ax.plot(te, roll_mean(yy, 60), lw=2.0, color=col)
    style_ax(ax)
    ax.set_title(cap, fontsize=12)
    ax.set_xlabel("Time (ns)" if a // 3 else "", fontsize=10.5)
    ax.set_ylabel(lab, fontsize=10.5)
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.text(0.985, 0.06, f"mean={yy.mean():.2f}", transform=ax.transAxes, ha="right",
            fontsize=9, color="#3C4852")
fig.suptitle("PprI(WT)·S1-ssDNA·Mn$^{2+}$   150 ns MD — 生产相能量与系综收敛 (pmemd.cuda mdout)",
             fontsize=14.5, fontweight="bold", y=0.965)
save(fig, "fig09_energy_convergence.png")

# ================================================================ FIG 10 RMSD 收敛
print("[fig10] RMSD convergence / 前半 vs 后半")
fig = plt.figure(figsize=(9.6, 5.6))
gs = fig.add_gridspec(2, 1, hspace=0.45, left=0.09, right=0.98, top=0.88, bottom=0.10)
ax = fig.add_subplot(gs[0])
half = NFR // 2
for key, col, lab in [("rmsd_prot", C_PROT, "Protein backbone"),
                      ("rmsd_dna", C_DNA, "ssDNA backbone")]:
    y = D[key]
    ax.plot(T, y, lw=0.5, color=col, alpha=0.25)
    ax.plot(T, roll_mean(y, 50), lw=2.2, color=col, label=lab)
    ax.axvline(T[half], ls="--", lw=0.9, color=GREY)
ax.text(T[half], ax.get_ylim()[1] * 0.98, " first half | second half ", ha="center",
        va="top", fontsize=8.5, color="#5A6670")
style_ax(ax); ax.set_ylabel("RMSD (Å)"); ax.set_xlabel("Time (ns)")
ax.legend(loc="upper left", ncol=2)
ax.set_title("时间窗平均 RMSD（收敛判据：后段不再漂移）", fontsize=12)
ax2 = fig.add_subplot(gs[1])
for key, col, lab in [("rmsd_prot", C_PROT, "protein"), ("rmsd_dna", C_DNA, "ssDNA")]:
    a = D[key][:half]; b = D[key][half:]
    for d, c, al in [(a, col, 0.34), (b, col, 1.0)]:
        kde = __import__("scipy.stats", fromlist=["gaussian_kde"]).gaussian_kde(d)
        xs = np.linspace(d.min() * 0.8, d.max() * 1.2, 220)
        ax2.plot(xs, kde(xs), color=c, lw=1.7, alpha=al)
        ax2.fill_between(xs, kde(xs), color=c, alpha=0.06 * (1 if al > 0.9 else 0.4))
    mu = f"{a.mean():.2f}→{b.mean():.2f} Å"
    ax2.text(0.98, 0.55 if key == "rmsd_prot" else 0.30, lab, transform=ax2.transAxes,
             ha="right", color=col, fontsize=10, fontweight="bold")
style_ax(ax2); ax2.set_xlabel("RMSD (Å)"); ax2.set_ylabel("Density")
ax2.set_title("前/后半程 RMSD 分布重叠度（虚线=后半，实线=前半；重合越好越收敛）", fontsize=12)
fig.suptitle("PprI(WT)·S1-ssDNA·Mn$^{2+}$   150 ns MD — 轨迹收敛性 (RMSD)", fontsize=14, y=0.965)
save(fig, "fig10_rmsd_convergence.png")

# ================================================================ FIG 11 RMSD-Rg 联合
print("[fig11] RMSD vs Rg joint")
fig = plt.figure(figsize=(7.6, 5.8))
ax = fig.add_subplot(111)
sc = ax.scatter(D["rmsd_prot"], D["rg_prot"], c=T, cmap="viridis", s=13,
                alpha=0.62, edgecolor="none", zorder=3)
cb = fig.colorbar(sc, ax=ax, pad=0.02); cb.set_label("Time (ns)", fontsize=11)
ax.set_xlabel("Protein backbone RMSD (Å)", fontsize=12.5)
ax.set_ylabel("Protein R$_g$ (Å)", fontsize=12.5)
style_ax(ax)
# 1σ 椭圆
from matplotlib.patches import Ellipse
mu = np.array([D["rmsd_prot"].mean(), D["rg_prot"].mean()])
cov = np.cov(D["rmsd_prot"], D["rg_prot"])
ev, evec = np.linalg.eigh(cov)
ang = np.degrees(np.arctan2(evec[1, -1], evec[0, -1]))
for n in (1, 2):
    w, h = 2 * n * np.sqrt(ev)
    ax.add_patch(Ellipse(mu, w, h, angle=ang, fill=False, ec=C_ACC, lw=1.8, ls="--", zorder=2))
stat_box(ax, f"⟨RMSD⟩ = {D['rmsd_prot'].mean():.2f} ± {D['rmsd_prot'].std():.2f} Å\n"
              f"⟨R$_g$⟩ = {D['rg_prot'].mean():.2f} ± {D['rg_prot'].std():.2f} Å", loc="upper left")
ax.set_title("构象采样：Rg–RMSD 联合分布（时间着色，虚线 = 1σ/2σ 椭圆）", fontsize=12.5)
save(fig, "fig11_rmsd_rg_joint.png")

# ================================================================ FIG 12 PCA 3D
print("[fig12] PCA 3D trajectory")
fig = plt.figure(figsize=(8.8, 6.6))
ax = fig.add_subplot(111, projection="3d")
cmap = plt.get_cmap("viridis")
c = [cmap(v) for v in np.linspace(0, 1, NFR)]
ax.scatter(D["pc1"], D["pc2"], D["pc3"], s=6, c=T, cmap="viridis", alpha=0.55, linewidth=0)
ax.plot(D["pc1"][::4], D["pc2"][::4], D["pc3"][::4], lw=0.6, color="#333", alpha=0.4)
# 首尾标记
ax.scatter([D["pc1"][0]], [D["pc2"][0]], [D["pc3"][0]], marker="o", s=90,
           color="#E8623C", edgecolor="white", zorder=6, label="start")
ax.scatter([D["pc1"][-1]], [D["pc2"][-1]], [D["pc3"][-1]], marker="*", s=220,
           color="#2B6CB0", edgecolor="white", zorder=6, label="end (150 ns)")
ax.set_xlabel("PC1 (48.0%)", fontsize=11); ax.set_ylabel("PC2 (13.4%)", fontsize=11)
ax.set_zlabel("PC3 (10.3%)", fontsize=11)
ax.view_init(elev=24, azim=-58)
ax.legend(loc="upper left", fontsize=10, frameon=False)
ax.set_title("主成分轨迹（PC1–PC3 三维投影，时间着色）", fontsize=12.5, pad=2)
# 手动 colorbar
sm = plt.cm.ScalarMappable(cmap="viridis", norm=mpl.colors.Normalize(0, 150))
cb = fig.colorbar(sm, ax=ax, shrink=0.55, pad=0.08); cb.set_label("Time (ns)", fontsize=10)
save(fig, "fig12_pca_3d.png")

# ================================================================ FIG 13 Ramachandran
print("[fig13] Ramachandran")
phi = np.loadtxt(os.path.join(RAW, "phi.dat"), skiprows=1)
psi = np.loadtxt(os.path.join(RAW, "psi.dat"), skiprows=1)
print(f"   phi shape {phi.shape} (frames x residues 2-254)")
# 全残基全帧 flatten (忽略极少数 NaN-like 大异常?直接裁剪)
ph = phi[:, 1:].ravel()          # 列0 -> 残基2
ps = psi[:, 1:].ravel()
ok = np.isfinite(ph) & np.isfinite(ps)
ph, ps = ph[ok], ps[ok]
fig = plt.figure(figsize=(8.2, 7.2))
ax = fig.add_subplot(111)
hb = ax.hexbin(ph, ps, gridsize=150, bins="log", cmap="Blues",
               extent=(-180, 180, -180, 180), mincnt=1)
cb = fig.colorbar(hb, ax=ax, pad=0.02, shrink=0.9); cb.set_label("log$_{10}$ count", fontsize=10.5)
# 关键残基散点
SEL = {67: "F88 (Core1 π)", 71: "H92 (HEXXH)", 149: "Y170 (Core1)", 232: "R253 (Core2)"}
colmap = {"F88 (Core1 π)": C_PROT, "H92 (HEXXH)": C_ACC, "Y170 (Core1)": C_PUR, "R253 (Core2)": C_TEAL}
for rseq, lab in SEL.items():
    col = phi[:, rseq - 2][::5]; ps2 = psi[:, rseq - 2][::5]
    ax.plot(col, ps2, lw=0.7, color=colmap[lab], alpha=0.55, zorder=4)
    ax.scatter(col[0], ps2[0], s=14, color=colmap[lab], zorder=5)
ax.set_xlim(-180, 180); ax.set_ylim(-180, 180)
ax.set_xticks(np.arange(-180, 181, 60)); ax.set_yticks(np.arange(-180, 181, 60))
ax.set_xlabel("φ (degrees)", fontsize=12.5); ax.set_ylabel("ψ (degrees)", fontsize=12.5)
ax.set_title("Ramachandran 分布（全部残基×全部帧）", fontsize=13)
handles = [plt.Line2D([], [], color=c, lw=2, label=l) for l, c in colmap.items()]
ax.legend(handles=handles, loc="lower right", fontsize=9, frameon=False,
          title="关键残基(seq)", title_fontsize=9)
save(fig, "fig13_ramachandran.png")

# ================================================================ FIG 14 DSSP
print("[fig14] DSSP secondary structure")
# dssp.dat: 每帧 1 行, 每残基 1 列 (code 整数). DSSP code 为标准约定:
#   0=coil/' ', 1=Extended(E), 2=Bridge(B), 3=3-10(G), 4=Alpha(H), 5=Pi(I), 6=Turn(T), 7=Bend(S)
dssp = np.loadtxt(os.path.join(RAW, "dssp.dat"), skiprows=1)[:, 1:].astype(int)
code2ss = {0: "coil", 1: "Extended", 2: "Bridge", 3: "3-10", 4: "Alpha", 5: "Pi",
           6: "Turn", 7: "Bend"}
print(f"   dssp {dssp.shape}; codes used: {sorted(np.unique(dssp).tolist())}")
# dssp_sum.dat: 第 0 列为 #Residue(1..254), 后 7 列 = Extended,Bridge,3-10,Alpha,Pi,Turn,Bend 占有率
sumo = np.loadtxt(os.path.join(RAW, "dssp_sum.dat"), skiprows=1)[:, 1:]
helix_codes = [3, 4, 5]      # G + H + I
sheet_codes = [1, 2]         # E + B
helix_frac = np.mean(np.isin(dssp, helix_codes), axis=1) * 100
sheet_frac = np.mean(np.isin(dssp, sheet_codes), axis=1) * 100
fig = plt.figure(figsize=(11.2, 6.8))
gs = fig.add_gridspec(2, 1, hspace=0.40, left=0.085, right=0.98, top=0.88, bottom=0.09,
                      height_ratios=[1, 0.72])
ax = fig.add_subplot(gs[0])
for yy, col, lab in [(helix_frac, C_PROT, "α-helix (H+G+I)"),
                     (sheet_frac, C_DNA, "β-strand (E+B)"),
                     (100 - helix_frac - sheet_frac, GREY, "coil/turn")]:
    ax.plot(T, yy, lw=1.6, color=col, label=lab)
style_ax(ax)
ax.set_ylim(-2, 102); ax.set_ylabel("Content (%)"); ax.set_xlabel("")
ax.legend(loc="lower left", ncol=3, fontsize=10)
ax.set_title("二级结构含量时程（每帧 DSSP）", fontsize=12.5)
ax.xaxis.set_major_locator(MultipleLocator(50))
ax2 = fig.add_subplot(gs[1])
# 每残基组成: sumo 列序 = Extended,Bridge,3-10,Alpha,Pi,Turn,Bend (idx 0-6)
hel = sumo[:, 2] + sumo[:, 3] + sumo[:, 4]   # 3-10 + Alpha + Pi
she = sumo[:, 0] + sumo[:, 1]                # Extended + Bridge
oth = 1 - hel - she
xr = np.arange(1, len(hel) + 1)
ax2.bar(xr, hel * 100, color=C_PROT, label="helix")
ax2.bar(xr, she * 100, bottom=hel * 100, color=C_DNA, label="strand")
ax2.bar(xr, oth * 100, bottom=(hel + she) * 100, color="#EDF0F3", label="coil")
ax2.set_xlim(0, 255); ax2.set_ylim(0, 100)
ax2.set_xlabel("Residue (seq)", fontsize=12); ax2.set_ylabel("Composition (%)")
for (a, b), col in [(HELIX_CAT, C_ACC), (LOOP_MOD, C_PUR)]:
    ax2.axvspan(a, b, color=col, alpha=0.14, lw=0)
ax2.text(HELIX_CAT[0] + 3, 92, "catalytic\nhelix", fontsize=8, color="#8A6D00")
ax2.text(LOOP_MOD[0] + 2, 92, "loop\n(modeled)", fontsize=8, color=C_PUR)
style_ax(ax2); ax2.legend(loc="lower left", ncol=3, fontsize=9)
ax2.set_title("每残基二级结构组成（DSSP 时间平均）", fontsize=12.5)
ax2.xaxis.set_major_locator(MultipleLocator(50))
fig.suptitle("PprI(WT)·S1-ssDNA·Mn$^{2+}$   150 ns MD — 二级结构 (DSSP)", fontsize=14, y=0.965)
save(fig, "fig14_dssp.png")

# ================================================================ FIG 15 界面接触数
print("[fig15] interface contacts")
nc = np.loadtxt(os.path.join(RAW, "nat_series.dat"), skiprows=1)
tn = nc[:, 0] * DT_NS * 5
native = nc[:, 1]
fig = plt.figure(figsize=(9.6, 4.4))
ax = fig.add_subplot(111)
ax.plot(tn, native, lw=0.7, color=C_PUR, alpha=0.4)
ax.plot(tn, roll_mean(native, 21), lw=2.3, color=C_PUR)
style_ax(ax)
ax.set_xlabel("Time (ns)", fontsize=12.5); ax.set_ylabel("# native contacts", fontsize=12.5)
ax.set_title("蛋白–DNA 界面接触数（重原子 ≤ 4.5 Å，相对首帧参考，每 0.5 ns）", fontsize=12.5)
ax.xaxis.set_major_locator(MultipleLocator(25))
stat_box(ax, f"mean = {native.mean():.0f} ± {native.std():.0f}\n"
             f"stable fraction = {native.mean()/native[0]*100:.1f}%", loc="upper right")
save(fig, "fig15_contacts_timeseries.png")

# ================================================================ FIG 16 接触热图
print("[fig16] contact map")
nr = np.loadtxt(os.path.join(RAW, "nat_res.dat"), skiprows=1)
prot_idx = nr[:, 0].astype(int)     # protein seq 1-254 (nativecontacts 只列出有接触的残基)
dna_idx = nr[:, 1].astype(int)      # DNA parm 残基号 255..278 (= nt1..24)
dna_min, dna_max = dna_idx.min(), dna_idx.max()
print(f"   protein res {prot_idx.min()}-{prot_idx.max()}, DNA res {dna_min}-{dna_max}")
DNA_PARM0 = 254                     # parm 编号 255 -> nt1 (即 nt = dna_idx - 254)
ntpos = dna_idx - DNA_PARM0
M = np.full((254, 24), np.nan)
M[prot_idx - 1, ntpos - 1] = nr[:, 2]      # TotalFrac (%)
missing = ~np.isin(np.arange(1, 25), ntpos)   # 从未接触的 nt -> 置 0 (游离前缀 TCATGAGC 等)
M[:, missing] = 0.0
fig = plt.figure(figsize=(10.4, 6.4))
ax = fig.add_subplot(111)
im = ax.imshow(M, aspect="auto", origin="lower", cmap="magma_r", vmin=0, vmax=100,
               extent=(0.5, 24.5, 0.5, 254.5))
ax.set_xlabel("ssDNA nucleotide (5′ → 3′)", fontsize=12.5)
ax.set_ylabel("Protein residue (seq)", fontsize=12.5)
ax.set_yticks(np.arange(1, 255, 20)); ax.set_xticks(np.arange(1, 25))
ax.set_xticklabels([f"{i}\n{S1[i-1]}" for i in range(1, 25)], fontsize=7.6)
# 关键核苷酸高亮
for npos in (10, 23):
    ax.axvline(npos, color="lime", lw=1.4, ls="--", alpha=0.8)
ax.text(10, 258, "dG10", color="lime", fontsize=9, ha="center", fontweight="bold")
ax.text(23, 258, "T23", color="lime", fontsize=9, ha="center", fontweight="bold")
# 关键蛋白残基行标
for rseq, (nm, grp) in KEY_RES.items():
    ax.plot(24.9, rseq, marker="<", ms=6, color=GROUP_COLOR[grp], clip_on=False)
cb = fig.colorbar(im, ax=ax, pad=0.015, shrink=0.9)
cb.set_label("Contact frequency (%)", fontsize=10.5)
cb.ax.tick_params(labelsize=9)
style_ax(ax, gridx=False, gridy=False)
ax.set_title("蛋白–DNA 残基接触图（距离 ≤ 4.5 Å，150 ns 平均；5′ 前缀 TCATGAGC 游离无接触）",
             fontsize=13, pad=26)
# 机制注释
top = nr[nr[:, 2].argsort()[::-1]][:3]
ax.text(0.0, -0.09, "▲ top 界面对：  " +
        "   ".join(f"{int(a)}·nt{int(b - DNA_PARM0)} ({c:.0f}%)" for a, b, c in
                  [(t[0], t[1], t[2]) for t in top]),
        transform=ax.transAxes, fontsize=9, color="#3C4852")
save(fig, "fig16_contact_map.png")

# ================================================================ FIG 17 界面氢键
print("[fig17] interface H-bonds")
rows = []
with open(os.path.join(RAW, "hb_all_avg.dat"), encoding="utf-8", errors="ignore") as fh:
    for ln in fh:
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.split()
        if len(f) < 7:
            continue
        acc, donh, don = f[0], f[1], f[2]
        m_acc = re.match(r"([A-Z0-9]+)_(\d+)@(.+)", acc)
        m_don = re.match(r"([A-Z0-9]+)_(\d+)@(.+)", don)
        if not m_acc or not m_don:
            continue
        ra, rd = int(m_acc.group(2)), int(m_don.group(2))
        is_if = (ra <= 254 < rd) or (rd <= 254 < ra)
        if not is_if:
            continue
        rows.append((float(f[4]), f[3], ra, rd, acc, don, f[5], f[6]))
rows.sort(reverse=True)
print(f"   {len(rows)} interfacial H-bonds; top: {rows[:3]}")
n_dna_donor = sum(1 for r in rows if r[3] > 254)
print(f"   DNA-donor count among all interface H-bonds: {n_dna_donor} / {len(rows)}")
top = rows[:15][::-1]
fig = plt.figure(figsize=(10.2, 5.6))
ax = fig.add_subplot(111)
names = []
fracs = []
for fr, nfr, ra, rd, acc, don, d, an in top:
    na = "DNA" if ra > 254 else acc.split("_")[0]
    nd = "DNA" if rd > 254 else don.split("_")[0]
    names.append(f"{don.split('_')[0]}{rd - 254 if rd > 254 else rd}·{acc.split('_')[0]}{ra - 254 if ra > 254 else ra}")
    fracs.append(fr * 100)
# occupancy 渐变着色: 0-100% → C_PROT 浅→深
import matplotlib as mpl2
cmap = mpl2.colors.LinearSegmentedColormap.from_list("occ_blue",
    ["#DCE6F2", C_PROT])
norm = mpl2.colors.Normalize(vmin=0, vmax=100)
cols = [cmap(norm(fr)) for fr in fracs]
bars = ax.barh(np.arange(len(top)), fracs, color=cols, alpha=0.95, height=0.62, ec="white", lw=0.4)
ax.set_yticks(np.arange(len(top)))
ax.set_yticklabels(names, fontsize=8.6)
ax.invert_yaxis()
ax.set_xlabel("Occupancy (%, 150 ns)", fontsize=12.5)
style_ax(ax, gridy=False)
for i, (fr, y) in enumerate(zip(fracs, np.arange(len(top)))):
    ax.text(fr + 1.2, y, f"{fr:.0f}%", va="center", fontsize=8.6, color="#3C4852")
ax.set_xlim(0, 108)
ax.set_title("界面氢键占用率 Top 15（蛋白供体→DNA 受体；色深=occupancy 高）", fontsize=13)
# colorbar 表示 occupancy 渐变
sm = mpl2.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = fig.colorbar(sm, ax=ax, pad=0.018, shrink=0.85)
cb.set_label("Occupancy (%)", fontsize=10)
cb.ax.tick_params(labelsize=8.5)
save(fig, "fig17_hbond_occupancy.png")

# ================================================================ FIG 18 口袋 RMSD
print("[fig18] catalytic pocket RMSD")
pk = read_dat(os.path.join(RAW, "pock_rmsd.dat"))
fig = plt.figure(figsize=(9.4, 4.4))
gs = fig.add_gridspec(1, 3, width_ratios=[5.2, 0.8, 0.8], wspace=0.06,
                      left=0.085, right=0.975, top=0.82, bottom=0.16)
ax = fig.add_subplot(gs[0, 0])
ax.plot(T, pk[:, 1], lw=0.9, color=C_PUR, alpha=0.42)
ax.plot(T, roll_mean(pk[:, 1], 50), lw=2.4, color=C_PUR)
style_ax(ax)
ax.set_xlabel("Time (ns)", fontsize=12.5); ax.set_ylabel("Pocket RMSD (Å)", fontsize=12.5)
ax.set_title("Mn$^{2+}$ 催化口袋（配位残基 H92/H96/E123 + Mn）", fontsize=12.5)
ax.xaxis.set_major_locator(MultipleLocator(50))
P.add_marginal(fig, gs[0, 1], pk[:, 1], C_PUR, xlabel="Density")
P.add_marginal(fig, gs[0, 2], pk[:, 1], C_PUR, xlabel="Density", orient="right")
stat_box(ax, f"mean = {pk[:,1].mean():.2f} ± {pk[:,1].std():.2f} Å", loc="upper right")
save(fig, "fig18_pocket_rmsd.png")

# ================================================================ FIG 19 DNA 端距
print("[fig19] DNA end-to-end")
ee = read_dat(os.path.join(RAW, "dna_ee.dat"))
fig = plt.figure(figsize=(9.4, 4.4))
gs = fig.add_gridspec(1, 3, width_ratios=[5.2, 0.8, 0.8], wspace=0.06,
                      left=0.085, right=0.975, top=0.82, bottom=0.16)
ax = fig.add_subplot(gs[0, 0])
ax.plot(T, ee[:, 1], lw=0.9, color=C_DNA, alpha=0.42)
ax.plot(T, roll_mean(ee[:, 1], 50), lw=2.4, color=C_DNA)
style_ax(ax)
ax.set_xlabel("Time (ns)", fontsize=12.5); ax.set_ylabel("End-to-end distance (Å)", fontsize=12.5)
ax.set_title("ssDNA 首尾距（5′C5′–3′O3′，柔性指标）", fontsize=12.5)
ax.xaxis.set_major_locator(MultipleLocator(50))
P.add_marginal(fig, gs[0, 1], ee[:, 1], C_DNA, xlabel="Density")
P.add_marginal(fig, gs[0, 2], ee[:, 1], C_DNA, xlabel="Density", orient="right")
stat_box(ax, f"mean = {ee[:,1].mean():.1f} ± {ee[:,1].std():.1f} Å", loc="upper right")
save(fig, "fig19_dna_end2end.png")

# ================================================================ FIG 20 RMSF 局部放大
print("[fig20] RMSF zoom (S167-Y217 region)")
fig = plt.figure(figsize=(10.4, 4.4))
ax = fig.add_subplot(111)
x0, x1 = 140, 205
xs = np.arange(1, len(D["rmsf_prot"]) + 1)
ax.plot(xs, D["rmsf_prot"], lw=1.6, color=C_PROT, alpha=0.35, zorder=2, label="full RMSF")
ax.plot(xs[x0:x1], D["rmsf_prot"][x0:x1], lw=2.4, color=C_PROT, zorder=3)
for (a, b), col, lab in [(HELIX_CAT, C_ACC, "catalytic helix (PDB 83-101)"),
                         (LOOP_MOD, C_PUR, "modeled loop (PDB 189-202)")]:
    if b < x0 or a > x1:
        continue
    ax.axvspan(max(a, x0), min(b, x1), color=col, alpha=0.14, lw=0)
    ax.text((a + b) / 2, 0.35, lab.split(" (")[0], color=col, fontsize=8,
            ha="center", va="bottom", fontweight="bold")
for rseq, (nm, grp) in KEY_RES.items():
    if x0 <= rseq <= x1:
        y0 = D["rmsf_prot"][rseq - 1]
        ax.scatter([rseq], [y0], s=52, color=GROUP_COLOR[grp], ec="white", zorder=6)
        ax.annotate(nm, (rseq, y0), (rseq, y0 + 1.0), ha="center", fontsize=9,
                    fontweight="bold", color=GROUP_COLOR[grp],
                    arrowprops=dict(arrowstyle="-", lw=0.9, color=GROUP_COLOR[grp]))
style_ax(ax)
ax.set_xlim(x0 - 2, x1 + 2)
ax.set_xlabel("Residue (seq)", fontsize=12.5); ax.set_ylabel("RMSF (Å)", fontsize=12.5)
ax.set_title("RMSF 局部放大：催化螺旋·建模环·双锁/核心区 (S167–Y217)", fontsize=12.5)
save(fig, "fig20_rmsf_zoom.png")

# ================================================================ FIG 21 核苷酸分解
print("[fig21] per-nucleotide MM-PBSA")
gb_rows = parse_decomp(os.path.join(RES, "FINAL_DECOMP_MMPBSA.dat")).get("GB", [])
nts = [r for r in gb_rows if r["loc"] == "L"]
nts.sort(key=lambda r: r["num"])
pos = np.array([r["num"] - 254 for r in nts])       # nt 1-24
tot = np.array([r["tot"] for r in nts])
fig = plt.figure(figsize=(11.4, 4.8))
ax = fig.add_subplot(111)
cols = [C_DNA if t < 0 else "#C9A2C8" for t in tot]
ax.bar(pos, tot, color=cols, alpha=0.9, width=0.78)
ax.axhline(0, lw=1.2, color=INK)
for i, t in enumerate(tot):
    ax.text(pos[i], t + (0.5 if t >= 0 else -0.5), f"{t:.1f}", ha="center",
            va="bottom" if t >= 0 else "top", fontsize=7.0, color="#3C4852")
ax.set_xticks(pos)
ax.set_xticklabels([f"{p}\n{S1[p-1]}" for p in pos], fontsize=8.6)
ax.set_xlabel("ssDNA nucleotide (5′ → 3′)", fontsize=12.5)
ax.set_ylabel("ΔG$_{nt}$ (kcal/mol, GB)", fontsize=12.5)
style_ax(ax)
ax.axvspan(9.5, 10.5, color="lime", alpha=0.16, lw=0)
ax.axvspan(22.5, 23.5, color="lime", alpha=0.16, lw=0)
ax.text(10, ax.get_ylim()[0] + (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.03,
        "dG10", color="#3E8E3E", fontsize=10, ha="center", fontweight="bold")
ax.text(23, ax.get_ylim()[0] + (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.03,
        "T23", color="#3E8E3E", fontsize=10, ha="center", fontweight="bold")
tot_nt = tot.sum()
ax.set_title(f"MM-PBSA per-nucleotide 分解（ΣΔG$_{{nt}}$ = {tot_nt:.1f} kcal/mol）", fontsize=13)
save(fig, "fig21_nucleotide_decomp.png")

# ================================================================ FIG 22 界面氢键二维谱
# 原 plan: 6 条关键界面氢键距离时序 (an3_hbdist.in + CHPC) -- 需要 CHPC 凭据
# 本地替身: 静态二维谱 occupancy vs avg-distance (来自 hb_all_avg.dat) +
#           注释关键键 (双锁/核心) + 角度着色 + 高亮 dG10/T23 关键键
print("[fig22] interface H-bond 2D spectrum (occupancy vs avg-distance)")
rows22 = []
with open(os.path.join(RAW, "hb_all_avg.dat"), encoding="utf-8", errors="ignore") as fh:
    for ln in fh:
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.split()
        if len(f) < 7:
            continue
        acc, donh, don = f[0], f[1], f[2]
        m_acc = re.match(r"([A-Z0-9]+)_(\d+)@(.+)", acc)
        m_don = re.match(r"([A-Z0-9]+)_(\d+)@(.+)", don)
        if not m_acc or not m_don:
            continue
        ra, rd = int(m_acc.group(2)), int(m_don.group(2))
        is_if = (ra <= 254 < rd) or (rd <= 254 < ra)
        if not is_if:
            continue
        rows22.append((float(f[4]), float(f[5]), float(f[6]),
                       ra, rd, acc, don, f[3]))
rows22.sort(reverse=True)
print(f"   {len(rows22)} interfacial H-bonds; top dist/angle summary ready")

# 关键键 (seq 蛋白 <-> nt): (donor_seq, acceptor_nt)  -> 标签
#   S167(seq146)-G17(nt17) ; N127(seq106)-G17(nt17) ; R207(seq186)-T15(nt15)
#   R253(seq232)-T23(nt23) ; S251(seq230)-T23(nt23) ; R220(seq199)-T22(nt22)
DNA_PARM0 = 254
KEY_BONDS = [
    (146, 271, "S167–G17", C_ACC),     # S167(seq146) Watson-Crick to G17(nt17,parm271)  90.7%
    (106, 271, "N127–G17", C_PROT),    # N127(seq106) Hoogsteen to G17                 76.0%
    (186, 269, "R207–T15", C_DNA),     # R207(seq186) to OP1 of T15(nt15,parm269)      87.0%
    (232, 278, "R253–T24", C_TEAL),    # R253(seq232) to OP1 of 3'T(nt24,parm278)      76.7%
    (230, 277, "S251–T23", C_PUR),     # S251(seq230) to O3' of T23(nt23,parm277)      74.0%
    (199, 276, "R220–T22", C_CPLX),    # R220(seq199) to OP1 of T22(nt22,parm276)      66.0%
]
KEY_LABELS = set((d, a) for d, a, *_ in KEY_BONDS)

fig = plt.figure(figsize=(10.6, 6.4))
gs = fig.add_gridspec(1, 1, left=0.085, right=0.985, top=0.90, bottom=0.11)
ax = fig.add_subplot(gs[0])
xs = np.array([r[1] for r in rows22])              # avg distance
ys = np.array([r[0] * 100 for r in rows22])        # occupancy %
ang = np.array([r[2] for r in rows22])             # avg angle (deg)
sc = ax.scatter(xs, ys, c=ang, cmap="viridis", s=18, alpha=0.7,
                vmin=120, vmax=180, edgecolor="none", zorder=3)
cb = fig.colorbar(sc, ax=ax, pad=0.018, shrink=0.92)
cb.set_label("Avg donor-H…acceptor angle (°)", fontsize=11)
ax.set_xlim(2.5, 4.0); ax.set_ylim(-3, 105)
ax.set_xlabel("Avg heavy-atom distance (Å)", fontsize=12.5)
ax.set_ylabel("H-bond occupancy (%, 150 ns)", fontsize=12.5)
# H-bond 几何阈值参考
ax.axvspan(2.5, 3.2, color=C_ACC, alpha=0.06, lw=0, zorder=1)
ax.axvspan(3.2, 3.5, color=C_ACC, alpha=0.03, lw=0, zorder=1)
ax.axvline(3.5, ls="--", lw=0.9, color=GREY, zorder=2)
ax.text(3.50, 100, "  3.5 Å cutoff (cpptraj)", ha="left", va="top",
        fontsize=9, color="#5A6670", zorder=2)
ax.text(2.85, 96, "strong", ha="center", va="top",
        fontsize=8.5, color="#3E8060", fontweight="bold")
ax.text(3.30, 96, "moderate", ha="center", va="top",
        fontsize=8.5, color="#3E8060", fontweight="bold")
# 高亮关键键
hl = []
for d_seq, a_parm, lab, col in KEY_BONDS:
    for r in rows22:
        if d_seq in (r[3], r[4]) and a_parm in (r[3], r[4]):
            hl.append((r, lab, col))
            break
for r, lab, col in hl:
    ax.scatter([r[1]], [r[0] * 100], s=140, facecolor="none",
               edgecolor=col, lw=2.0, zorder=4)
    ax.annotate(lab, (r[1], r[0] * 100),
                xytext=(r[1] + 0.04, r[0] * 100 - 4.5),
                fontsize=9, color=col, fontweight="bold",
                arrowprops=dict(arrowstyle="-", lw=0.8, color=col))
# 整体标注
ax.set_title("界面 H 键二维谱：occupancy × 距离 × 角度（关键双锁/核心键高亮）",
             fontsize=13, pad=6)
# 上右框去除
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.tick_params(direction="out")
ax.grid(True, lw=0.3, color=GRIDC, alpha=0.5)
stat_box(ax, f"n = {len(rows22)} interface H-bonds\n"
             f"median occ. = {np.median(ys):.0f} %\n"
             f"strong (<3.2 Å) = {(xs < 3.2).sum()}",
         loc="upper right")
fig.suptitle("PprI(WT)·S1-ssDNA·Mn$^{2+}$   150 ns MD — 界面 H 键网络静态谱", fontsize=14, y=0.972)
save(fig, "fig22_hbond_timeseries.png")
print(f"   [fig22] fallback used (static 2D spectrum); an3-based timeseries needs CHPC")

print("[done] fig09-22 all saved to", OUT)
