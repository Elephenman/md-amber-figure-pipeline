#!/usr/bin/env python3
"""生成 22 图集 HTML 画廊 — figures_pretty/index.html (深色科研风单文件)"""
import os, html

ITEMS = [
    # (fig_num, fname, group, title, key_value, caption)
    ("01","fig01_rmsd.png","基础","蛋白 + ssDNA RMSD 时程",
     "prot 1.93±0.36 Å  ·  ssDNA 7.41±1.41 Å",
     "骨架 RMSD 1500 帧 × 100 ps；浅色=原曲线，深色=50 帧滑动平均。蛋白骨架已收敛（后段不再漂移），ssDNA 较高反映 5′ 延伸 9 nt 自由摆动。"),
    ("02","fig02_rmsf.png","基础","每残基 RMSF（蛋白 254 aa）",
     "⟨RMSF⟩ 1.06±0.74 Å",
     "B-factor 等价 CA RMSF。催化螺旋 (seq 62-80) 与建模环 (168-181) 柔性高，与已知功能区一致。"),
    ("03","fig03_rg.png","基础","蛋白回转半径 Rg",
     "⟨Rg⟩ 19.94±0.23 Å",
     "150 ns 内 Rg 标准差仅 0.23 Å → 整体拓扑无膨胀/塌缩。"),
    ("04","fig04_sasa.png","基础","SASA 复合物 vs 分离",
     "buried 3445±264 Å²",
     "5 状态时程：cplx_total, prot_in_cplx, prot_iso, dna_in_cplx, dna_iso；buried = (prot_iso - prot_in_cplx) + (dna_iso - dna_in_cplx)。<em>关键坑</em>：复合物 SASA 不能用 prot_iso + dna_iso - cplx（重叠抵消会让面积虚增）。"),
    ("05","fig05_pca.png","集合","PC1-PC2 + 自由能形貌 FEL",
     "PC1 48% / PC2 13%",
     "CA 协方差前 2 主成分 + 100×100 bin KDE FEL。色=自由能 kcal/mol，单能谷 ΔG ≈ 13 kcal/mol 高度。"),
    ("06","fig06_fel.png","集合","FEL 等高线 + 极小点",
     "1 个主要能谷",
     "聚焦能谷在 PC1≈0 / PC2≈0，轨迹主要在此盆地内采样。"),
    ("07","fig07_mmpbsa.png","能量","MM-PBSA 逐残基 ΔG 分解",
     "top 贡献在双锁/核心区",
     "GB 模型（不显式水），ΔG = EEL + VDW + ΔG_sol。关键残基贡献明显（深蓝/深红），与机制图谱对应。"),
    ("08","fig08_mn_coordination.png","能量","Mn²⁺ 三齿配位距离",
     "H71 2.23±0.08 / H75 2.25±0.08 / E102a 2.77±0.74 Å",
     "<strong>HEXXH-E 三齿配位保持稳定</strong>（H92/H96 严格保持，E123 较软）。Mn²⁺ 在 150 ns 内未解离，配位几何与晶体一致。"),
    ("09","fig09_energy_convergence.png","收敛","能量 + 温度 + 密度收敛",
     "T=303 K  ·  P~1 bar  ·  ρ~1.0 g/cm³",
     "TEMP(K)/PRESS(bar)/Density/Etot/Ep/Ek 共 6 子图，60 帧滑动平均。NPT 系综稳定收敛。"),
    ("10","fig10_rmsd_convergence.png","收敛","RMSD 前半 vs 后半 KDE",
     "KDE 重叠好",
     "前 750 帧 vs 后 750 帧的 RMSD KDE 叠加，曲线几乎重合 → 轨迹已充分平衡，可作统计采样。"),
    ("11","fig11_rmsd_rg_joint.png","收敛","RMSD-Rg 联合散点 (时间着色)",
     "2σ 椭圆内主要密度",
     "RMSD×Rg 相空间散点，色=时间(ns)。1σ/2σ 协方差椭圆，绝大部分轨迹落在椭圆内。"),
    ("12","fig12_pca_3d.png","收敛","PC1-PC3 三维轨迹",
     "首尾 start / end (150 ns)",
     "3D 投影，色=时间。低维投影保留大部分构象动力学。"),
    ("13","fig13_ramachandran.png","构象","Ramachandran 全部残基×全部帧",
     "log10 count 密度",
     "1500 帧 × 253 个 φψ 对 hexbin；关键残基 F88 / H92 / Y170 / R253 轨迹叠加（高亮彩线）。"),
    ("14","fig14_dssp.png","构象","DSSP 二级结构（时程+残基组成）",
     "⟨α-helix⟩ 42%  ·  ⟨β-strand⟩ 22%  ·  ⟨coil⟩ 35%",
     "<em>修复点</em>：硬编码 DSSP code 映射 + <code>sumo[:,1:]</code> 切列。催化螺旋 (res 62-80) 与建模环 (168-181) 在残基组成图中清晰可见。"),
    ("15","fig15_contacts_timeseries.png","界面","蛋白-DNA 界面接触数时程",
     "⟨#contacts⟩ 634±108 (32% of initial)",
     "nativecontacts ≤4.5 Å per residue，初始 1982 → 末期 ~634（20 帧滑动平均）。32% 保留率反映核心区持续接触、外围接触演化。"),
    ("16","fig16_contact_map.png","界面","残基-核苷酸接触热图",
     "★ S1 5′ TCATGAGC 9 nt 完全无接触",
     "<em>修复点</em>：固定 <code>DNA_PARM0=254</code> 修正 nt 偏移。<strong>关键发现</strong>：5′ 端 9 nt (TCATGAGC) 从 MD 起始到 150 ns 始终无接触 = 游离尾巴；接触集中 nt10-24（G17 关键锚点 / T23 端点）。"),
    ("17","fig17_hbond_occupancy.png","界面","界面氢键占用率 Top 15 (蛋白→DNA)",
     "S167-G17 90.7%  ·  R207-T15 87.0%  ·  R253-T24 76.7%",
     "<em>修复点</em>：top-15 全部为蛋白供体（0 DNA 供体），改 occupancy 渐变着色。前 3 涵盖双锁/核心 1/核心 2 三组关键键。"),
    ("18","fig18_pocket_rmsd.png","局部","Mn 催化口袋 RMSD",
     "⟨pocket RMSD⟩ 0.54±0.08 Å",
     "残基 (H92/H96/E123 + Mn) 配位集合 RMSD，150 ns 内极稳定（亚 Å），支持\"催化口袋刚性\"假设。"),
    ("19","fig19_dna_end2end.png","局部","ssDNA 首尾距 (5′C5′–3′O3′)",
     "⟨EE⟩ 47.6±6.2 Å",
     "24 nt 链端到端距离，每 0.1 ns 一次。C5′ 起算因 5′DT5 无 P 原子。"),
    ("20","fig20_rmsf_zoom.png","局部","RMSF 局部放大 (S167-Y217 区)",
     "max RMSF 5.44 Å @ res 174",
     "seq 140-205 局部 RMSF，标注 KEY_RES（双锁/核心区/建模环）。res 174 (PDB 195) 在建模环内 = 最柔性残基。"),
    ("21","fig21_nucleotide_decomp.png","机制","MM-PBSA 逐核苷酸分解",
     "G17 / T23 双极小",
     "S1 24 nt 逐位 ΔG，dG10 (pos10) 与 T23 (pos23) 高亮。负值=有利接触。"),
    ("22","fig22_hbond_timeseries.png","机制","界面 H 键 动态时序 (6 关键键 150 ns)",
     "核心1 全程稳定  ·  核心2 ~120 ns 后解离",
     "6 条关键界面 H 键逐帧 (0.1 ns) 距离 + 成键栅栏 (CHPC cpptraj V6.4.4 an3_hbdist.in)。<br><b>科学故事</b>：核心 1 (S167–G17 / N127–G17) + Patch1 (R207–T15) 全程 71–91% 维持；核心 2 (R253–T24 / S251–T23 / R220–T22) 前 75 ns 尚未建立，<b>~120 ns 后 R220/R253/S251 同时抬升 &gt; 4 Å</b>，栅栏灰条=解离/展开。"),
    ("23","fig23_dccm.png","机制","蛋白运动互相关矩阵 DCCM (254×254)",
     "vlim ±0.79  ·  cat×Cterm ⟨C⟩=−0.11",
     "蛋白 CA 原子 rms-fit 后两两互相关 (CHPC <code>matrix correl</code>)。对角 = 自相关。HEXXH 催化螺旋 (62-80) 与 C 端核酸结合区 (180-254) 块内强正相关；催化 × C 端 = −0.11 (微弱反相关) → 提示跨域偶联。"),
    ("24","fig24_mmpbsa_prot_residue.png","能量","MM-PBSA 蛋白 per-residue ΔG 分解",
     "ΣΔG = −126.2 kcal/mol  ·  R207 = −13.0 top1",
     "GB 范式下 254 蛋白残基的 vdW+EEL+EGB+ESURF 贡献。深蓝 = 稳定结合；浅红 = 不利。GAF 域/蛋白酶域背景淡化。稳定 top: <b>R207 −13.0</b>（Patch1 主锚）, 168 −10.2, R250 −9.4, R220 −9.0, R253 −8.8。"),
    ("25","fig25_hbond_occupancy_dist.png","界面","界面 H 键 occupancy 分布",
     "n=372  ·  ≥60% 仅 7 键  ·  长尾分布",
     "界面 (蛋白×DNA) H 键 occupancy 直方图 + 累积曲线（左）+ occ×dist 密度散射（右）。“少数键承担主要结合”——≥60% 持久键 7 个，&lt;20% 瞬时键 336。"),
    ("26","fig26_contacts_hbond_dual.png","界面","界面接触 + H 键 双通道时序",
     "native 634±108  ·  H-bond 269±9",
     "nativecontacts (native + non-native) 与 hbond 总数 (HBALL[UU]) 双 y 轴叠加，7.5 ns 滚动均值。stride 5 = 300 帧。两条曲线 <b>同步波动</b> 证明接触网络与 H 键网络是同一物理事件的两个表现。"),
    ("27","fig27_dssp_timeline.png","构象","DSSP 二级结构时间演化栅栏 (1500 帧 × 254 残基)",
     "α 占比 0.41±0.01  ·  loop 168-181 高度动态",
     "DSSP 逐帧 (150 ns) 颜色编码二级结构 (红 α / 黄 β / 绿 turn / 蓝 3-10 / 青 bend / 紫 π / 灰 coil)。催化螺旋 (62-80) 主体 α；无序环 (168-181) 高度混合。"),
    ("28","fig28_rmsf_vs_bfactor.png","局部","MD RMSF vs 晶体 B-factor 校验",
     "Pearson r = 0.62 (n=240)",
     "晶体 CA B-factor 经 <code>√(3B/8π²)</code> 当量化与 MD RMSF 同轴对比。无序环 (168-181) 在晶体中缺失，MD 中表现为 5.4 Å 最高峰 (建模补出)。Pearson r=0.62 验证 MD 柔性场与晶体 B 良好一致。"),
    ("29","fig29_dna_morphology.png","构象","ssDNA 形态三指标 (端距 / Rg / 弯曲)",
     "EE 47.6±6.2  ·  Rg 21.4±1.1  ·  bend 85±15°",
     "三面板 (3 ns 滚动均值): a 端到端距离; b 回旋半径; c 弯曲角 (5′-中-3′ C3′)。<b>~120 ns 后三指标同步下降</b> (EE→37 Å, Rg→19.5 Å, bend→60°): ssDNA 在界面向更紧凑、更弯曲的构象重排——与 fig22 关键 H 键解离时间精准对应。"),
    ("30","fig30_mn_rdf.png","机制","Mn²⁺ 配位环境径向分布 (RDF)",
     "第一壳层 1.95 Å (g=6.53)  ·  第二壳层 3.15 Å",
     "蛋白原子围绕 Mn²⁺ 的 g(r)。第一壳层 1.95 Å 强峰 = Mn 与配位原子 (H71/H75/E102 的 N/O) 直接配位；第二壳层 3.15 Å = 第二层配位残基。参考线 H71 (2.23), H75 (2.25), E102 (2.77) 配位距离。"),
]

GROUPS = ["基础","集合","能量","收敛","构象","界面","局部","机制"]
GROUP_COLORS = {
    "基础":"#2B6CB0", "集合":"#7E57C2", "能量":"#10A37F",
    "收敛":"#3A8FB7", "构象":"#8E44AD", "界面":"#E8623C",
    "局部":"#3A8FB7", "机制":"#D69E2E",
}

def esc(s): return html.escape(str(s))

cards = []
for fig, fn, grp, title, kv, cap in ITEMS:
    color = GROUP_COLORS.get(grp, "#2B6CB0")
    cards.append(f"""
    <a class="card" id="fig{fig}" href="{fn}" target="_blank">
      <div class="card-head">
        <span class="fig-no" style="background:{color}">fig{fig}</span>
        <span class="fig-grp" style="color:{color}">{esc(grp)}</span>
      </div>
      <div class="thumb"><img src="{fn}" loading="lazy" alt="fig{fig} {esc(title)}"></div>
      <div class="title">{esc(title)}</div>
      <div class="kv" style="border-left:3px solid {color}">{esc(kv)}</div>
      <div class="cap">{cap}</div>
    </a>""")

toc = " ".join(
    f'<a class="toc-grp" style="color:{GROUP_COLORS[g]}" href="#grp-{i}">{g}</a>'
    for i, g in enumerate(GROUPS) if any(grp==g for _,_,grp,_,_,_ in ITEMS)
)

grouped_cards = ""
for gi, g in enumerate(GROUPS):
    rows = [c for c,(_,_,grp,_,_,_) in zip(cards, ITEMS) if grp==g]
    if not rows: continue
    grouped_cards += f'<h2 id="grp-{gi}" class="grp-title" style="border-bottom:2px solid {GROUP_COLORS[g]}">{esc(g)} <span class="grp-count">{len(rows)} 张</span></h2><div class="grid">{"".join(rows)}</div>'

html_doc = f"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<title>PprI(WT)·S1-ssDNA·Mn²⁺ 150 ns MD 图集 (22 张)</title>
<style>
  :root {{
    --bg: #0F1419;  --panel: #161D27;  --border: #1E2837;
    --ink: #E2E8F0;  --dim: #8A95A5;  --accent: #2B6CB0;
  }}
  * {{ box-sizing: border-box; }}
  body {{ margin: 0; background: var(--bg); color: var(--ink);
          font: 14px/1.55 -apple-system,"Microsoft YaHei",sans-serif; }}
  header {{ padding: 36px 48px 24px; background: linear-gradient(180deg, #1A2230, #0F1419);
            border-bottom: 1px solid var(--border); }}
  header h1 {{ margin: 0 0 8px; font-size: 26px; color: #FFF; font-weight: 600; letter-spacing: 0.5px; }}
  header h1 small {{ color: #E8623C; font-weight: 400; margin-left: 8px; font-size: 18px; }}
  .meta {{ color: var(--dim); font-size: 13px; line-height: 1.8; }}
  .meta code {{ background: #1E2837; color: #E8623C; padding: 1px 6px; border-radius: 3px; font-size: 12px; }}
  .meta a {{ color: #7AB7E8; text-decoration: none; }}
  .meta a:hover {{ text-decoration: underline; }}
  .toc {{ padding: 14px 48px; background: #131A24; border-bottom: 1px solid var(--border);
         position: sticky; top: 0; z-index: 9; }}
  .toc-grp {{ display: inline-block; margin-right: 18px; padding: 4px 10px; font-size: 13px;
              font-weight: 600; text-decoration: none; border-radius: 3px;
              background: #1E2837; transition: background .15s; }}
  .toc-grp:hover {{ background: #2A3648; }}
  main {{ padding: 24px 48px 60px; max-width: 1700px; margin: 0 auto; }}
  .grp-title {{ margin: 36px 0 16px; padding-bottom: 8px; font-size: 20px; font-weight: 600; }}
  .grp-count {{ font-size: 13px; color: var(--dim); font-weight: 400; margin-left: 6px; }}
  .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(440px, 1fr));
           gap: 18px; }}
  .card {{ display: block; background: var(--panel); border: 1px solid var(--border);
           border-radius: 6px; padding: 14px 16px 18px; text-decoration: none;
           color: var(--ink); transition: transform .15s, border-color .15s, box-shadow .15s; }}
  .card:hover {{ transform: translateY(-2px); border-color: #3A4960;
                  box-shadow: 0 6px 18px rgba(0,0,0,0.4); }}
  .card-head {{ display: flex; align-items: center; gap: 10px; margin-bottom: 8px; }}
  .fig-no {{ display: inline-block; padding: 2px 8px; border-radius: 3px;
             color: #FFF; font-family: ui-monospace,"SF Mono",Consolas,monospace;
             font-size: 12px; font-weight: 600; letter-spacing: 0.5px; }}
  .fig-grp {{ font-size: 12px; font-weight: 600; letter-spacing: 1px; }}
  .thumb {{ display: flex; justify-content: center; background: #FFFFFF;
            border-radius: 4px; margin: 8px 0 10px; min-height: 180px; align-items: center; }}
  .thumb img {{ max-width: 100%; max-height: 320px; display: block; }}
  .title {{ font-size: 14.5px; font-weight: 600; margin: 6px 0 8px; color: #F7FAFC; }}
  .kv {{ background: #0E141C; padding: 6px 10px; margin: 0 0 10px; font-size: 12.5px;
         color: #C4D0DE; font-family: ui-monospace,"SF Mono",Consolas,monospace; border-radius: 0 3px 3px 0; }}
  .cap {{ font-size: 12.8px; line-height: 1.65; color: #B6BFCB; }}
  .cap strong {{ color: #E8623C; }}
  .cap code {{ background: #1E2837; color: #C4D0DE; padding: 1px 5px; border-radius: 2px; font-size: 11.8px; }}
  .cap em {{ color: #B0B7C2; font-style: normal; font-weight: 600; }}
  footer {{ padding: 30px 48px; color: var(--dim); font-size: 12px; text-align: center;
            border-top: 1px solid var(--border); }}
  footer a {{ color: #7AB7E8; text-decoration: none; }}
</style>
</head>
<body>
<header>
  <h1>PprI(WT)·S1-ssDNA·Mn²⁺ 150 ns MD 图集<small> · 22 张</small></h1>
  <div class="meta">
    数据源：CHPC 4090 真实 MD（prod.nc 1500 帧 × 100 ps）<br>
    作图脚本：<code>plot_pretty_figs.py</code> + <code>plot_pretty_figs2.py</code>（22 张推文级 matplotlib 渲染）<br>
    工具链：AmberTools 24 cpptraj · MMPBSA.py (GB) · numpy/scipy/matplotlib<br>
    风格：<a href="https://github.com/Elephenman/md-amber-figure-pipeline" target="_blank">md-amber-figure-pipeline</a>（public repo）
  </div>
</header>
<div class="toc">{toc}</div>
<main>
{grouped_cards}
</main>
<footer>
  编号约定：蛋白 seq 1-254（parm = seq + 21），DNA parm 255-278 = nt 1-24；S1 = TCATGAGCAGTTTTTTGTTTTTTT (24 nt)
  · 关键发现：5′ TCATGAGC 9 nt 完全无接触（= 游离尾巴）；S167-G17 H 键 90.7% 持续 = 双锁锚点
</footer>
</body>
</html>"""

out_path = r"A:\Data\设计蛋白\PprI_ssDNA_design\md_specificity\results\WT__S1\figures_pretty\index.html"
with open(out_path, "w", encoding="utf-8") as fh:
    fh.write(html_doc)
print(f"[written] {out_path}  ({len(html_doc):,} bytes)")
