# md-amber-figure-pipeline

端到端 **Amber MD 后分析 + 推文级美化作图** 流水线。从 cpptraj/MMPBSA 跑出的 ASCII 数据出发，渲染 22 张高质量图（150 ns 蛋白·核酸 MD 标配）。

> 模板来自 **PprI·ssDNA·Mn²⁺ 150 ns MD** 项目验证（[Elephenman/PprI_ssDNA_design](https://github.com/Elephenman/PprI_ssDNA_design)）。本仓库只装脚本与示例输出，不装轨迹/原始数据。

## 22 张图一览

| 类别 | 图号 | 内容 |
|---|---|---|
| 基础 | fig01-04 | RMSD / RMSF / Rg / SASA (含 cplx vs isolated → buried 面积) |
| 集合 | fig05-06 | PC1-2 + FEL |
| 能量 | fig07-08 | MM-PBSA 逐残基 + 金属三齿配位距离 |
| 收敛 | fig09-12 | 能量时程 + RMSD 收敛 + Rg-RMSD 联合 + PC1-3 三维 |
| 构象 | fig13-14 | Ramachandran + DSSP 二级结构 |
| 界面 | fig15-17 | 接触数时程 + 接触热图 + 氢键占用率 Top 15 |
| 局部 | fig18-20 | 催化口袋 RMSD + DNA 首尾距 + RMSF 局部放大 |
| 机制 | fig21-22 | MM-PBSA 逐核苷酸 + 6 条关键界面氢键距离时序 |

完整图集见 `examples/fig01..fig22.png`。

## 数据布局（必须）

```
<your-system>/
  plot_pretty_figs.py
  plot_pretty_figs2.py
  cpptraj/
    an1_struct.in
    an2b_contacts.in
    an2c_pock.in
    an2d_ee.in
    an3_hbdist.in
  results/<your-system>/
    prod.out                       # Amber mdout (fig09)
    rmsd_prot.dat, rmsd_dna.dat    # cpptraj rms
    rg_prot.dat, rg_dna.dat
    rmsf_prot.dat, rmsf_dna.dat
    sasa_cplx.dat, sasa_prot.dat, sasa_dna.dat
    sasa_prot_iso.dat, sasa_dna_iso.dat
    pca_all.dat, pca_proj.dat
    mn_H*.dat, mn_E*.dat            # 金属配位 (如适用)
    mn_rdf.dat
    mmpbsa.log
    FINAL_DECOMP_MMPBSA.dat        # MM-PBSA 逐残基分解
    cpptraj_raw/                   # 6 个 an*.in 输出
      dssp.dat, dssp_sum.dat
      phi.dat, psi.dat
      nat_res.dat, nat_series.dat
      hb_all_avg.dat               # (可选, fig17 用)
      pock_rmsd.dat, dna_ee.dat
      hbd_*.dat                    # an3 输出 6 个距离文件
  figures_pretty/                  # 输出 PNG 目录
```

## 上游数据获取

### 1. cpptraj 原始数据采集
```bash
# 假设在 CHPC 集群登录节点 (10.202.94.52:20009) 或本地 AmberTools 环境
cd ~/md/out/<your-system>/
for f in an1_struct an2b_contacts an2c_pock an2d_ee an3_hbdist; do
  cpptraj -i ~/md/cpptraj/${f}.in > ${f}.log 2>&1
done
```
**修改 `cpptraj/an*.in` 里的 parm/trajin 路径到你的实际产物**。

### 2. MM-PBSA 逐残基分解（如已有可跳过）
```bash
MMPBSA.py -O -i mmpbsa.in -o FINAL_RESULTS.dat \
  -eo FINAL_DECOMP_MMPBSA.dat -sp complex.parm7 -cp complex.parm7 \
  -rp protein.parm7 -lp ligand.parm7 -y prod.nc
```

## 一键作图

```bash
# 1) 修改 plot_pretty_figs.py 顶部 6 个常量
#    CASE, RES, OUT, D, NFR, DT_NS, T
#    改 KEY_RES 残基映射 (seq -> PDB 名, group)
#    改 S1 = "<24nt 序列>" (仅核酸体系)

# 2) 跑基础 8 张 (fig01-08)
python scripts/plot_pretty_figs.py

# 3) 跑扩展 14 张 (fig09-22)
python scripts/plot_pretty_figs2.py
```

输出自动到 `figures_pretty/fig01..fig22.png`，1500 帧/150ns 约 25 秒。

## 关键设计约定

| 项目 | 约定 |
|---|---|
| 编号 | 蛋白 seq 1..N（parm 编号 = seq 编号 + 21），DNA 起始 255 = nt1 |
| 体系规模 | 30-100k 原子 |
| 帧 | 1000-3000 帧（默认 1500 × 0.1 ns = 150 ns） |
| 字体 | Microsoft YaHei + DejaVu Sans（已 fallback 硬编码） |
| 配色 | C_PROT=#2B6CB0, C_DNA=#E8623C, C_ACC=#10A37F, C_PUR=#8E44AD, C_TEAL=#3A8FB7 |

## 已知坑（必读）

- **DSSP code**：硬编码 `{0:coil, 1:E, 2:B, 3:G, 4:H, 5:I, 6:T, 7:S}`；`dssp_sum.dat` 第 0 列是 `#Residue`，**先 `[:,1:]`**
- **nat_res.dat DNA 偏移**：固定 254（parm 255=nt1），不要用 `dna_min-1`
- **hb_all_avg 界面过滤**：`is_if = (ra<=254<rd) or (rd<=254<ra)`
- **Microsoft YaHei 缺 `⟨⟩`**：用 `<>` 或 `mean()` 文本
- **`hbond` 双 mask** 必须用 `,` 合并 `:1-254,:255-278`，否则报错；但合并后 series 只有 `HBALL[UU]` 总数 → 逐键距离用 `an3_hbdist.in`
- **孤立 SASA**：必须 `strip + parmstrip` 同步剥除，否则=0
- **Boltz vs AF3 ipTM 系统偏差** 0.05-0.1，**同平台比较**
- **跨预测 RMSD** 必须先 Kabsch 叠加

详见 `scripts/plot_pretty_figs.py` 顶部注释与每个 fig 的 `print()` 输出。

## 依赖

- Python 3.10+
- numpy ≥ 1.24
- scipy ≥ 1.10
- matplotlib ≥ 3.7
- AmberTools 24+ (跑 cpptraj / MMPBSA.py)

`pip install numpy scipy matplotlib` 即可，不需要 GPU。

## 复现示例

仓库 `examples/` 下 22 张 PNG 是 **PprI·ssDNA·Mn²⁺ WT 体系** 150 ns 真实 MD（蛋白 seq1-254, ssDNA nt1-24, S1 序列 `TCATGAGCAGTTTTTTGTTTTTTT`）的输出，可作图质量参考。

## 引用

```
@software{ele_2026_md_pipeline,
  author = {Ye, Yongfeng (Elephenman)},
  title  = {md-amber-figure-pipeline: Amber MD 后分析 + 推文级美化作图 22 张图集},
  year   = {2026},
  url    = {https://github.com/Elephenman/md-amber-figure-pipeline}
}
```

## License

MIT
