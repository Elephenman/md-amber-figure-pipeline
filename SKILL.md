---
name: md-amber-figure-pipeline
description: 端到端的 Amber MD 后分析 + 推文级美化作图流水线。从 cpptraj/MMPBSA 跑出的 ASCII 数据出发，渲染 30 张高质量图（150 ns 蛋白·核酸 MD 标配）。适用于蛋白-RNA/DNA/小分子复合物 MD 结果整理。触发词：MD 漂亮图、MD 美化作图、Amber 轨迹分析、cpptraj 出图、推文级 MD 图、md 30 图、WT MD 出图、150 ns MD 作图。
agent_created: true
---

# MD Amber 美化作图流水线（PprI·ssDNA 150 ns 模板）

## 一句话
给一份 Amber 跑完的复合物 MD（parm7 + prod.nc + MMPBSA 结果），按 8 个 cpptraj 脚本采集分析数据，喂给 3 个 Python 脚本，**自动产出 30 张** 推文级 PNG。

## 适用对象
- 蛋白·核酸 / 蛋白·小分子 / 蛋白·蛋白复合物 MD
- 50-300 ns 段，1000-3000 帧（默认 1500 帧 = 150 ns × 0.1 ns/帧）
- 体系规模 30-100k 原子
- 残基编号约定：seq 1..N（parm 编号 = seq 编号 + 21，DNA 起始 255 = nt 1）
- **前提**：有真 MD 产物（轨迹 parm7 + prod.nc + MMPBSA `FINAL_DECOMP_MMPBSA.dat`）

## 文件结构
```
<system>/
  md_specificity/
    plot_pretty_figs.py          # 图集 v1 (fig01-08) — RMSD/RMSF/Rg/SASA/PCA/FEL/MM-PBSA/Mn
    plot_pretty_figs2.py         # 图集 v2 (fig09-22) — 能量收敛/二级结构/界面接触/氢键/pocket/DNA端距/MM-PBSA 逐核苷酸
    plot_pretty_figs3.py         # 图集 v3 (fig23-30) — DCCM/per-residue ΔG/H键occ/双通道/DSSP timeline/RMSF vs B/DNA形态/Mn RDF
    cpptraj/
      an1_struct.in              # dssp + multidihedral (phi/psi)
      an2_interface.in           # hbond (会报错,保留作记录)
      an2b_contacts.in           # nativecontacts + hbond series  (主接触分析)
      an2c_pock.in               # 催化口袋 RMSD
      an2d_ee.in                 # DNA 首尾距
      an3_hbdist.in              # 6 条关键界面氢键距离时序 (fig22 完整版数据源)
      an4_dynamics.in            # DNA 弯曲角 dna_bend (fig29 末面板)
      an4b_dccm_fit.in           # DCCM (fig23) — 必须先 rms fit 去刚体
      an5_rdf.in                 # Mn²⁺ 配位 RDF (fig30) — mask1=蛋白, mask2=MN
    results/<system>/
      prod.out                   # Amber mdout, fig09 能量收敛用
      *.dat                      # cpptraj 原始输出
      FINAL_DECOMP_MMPBSA.dat    # MM-PBSA 逐残基分解
      bfactor_prot.dat           # 8SLN 晶体链A CA B-factor (fig28) — 本地从 RCSB 8SLN.pdb 抽
      cpptraj_raw/               # 8 个 an*.in 输出归集
        dssp.dat, dssp_sum.dat, phi.dat, psi.dat
        nat_res.dat, nat_series.dat
        hb_all_avg.dat           # (图17/25 界面氢键占用率来源)
        pock_rmsd.dat, dna_ee.dat
        hbd_*.dat                # (6 条 an3 输出, fig22 完整时序版必备)
        dna_bend.dat             # (an4 输出, fig29 c 面板)
        dccm_prot.dat            # (an4b 输出, fig23)
      figures_pretty/fig01..fig30.png
```

## 标准流程（5 步）
1. **跑 cpptraj 采集原始数据**（在 CHPC 4090 集群登录节点或本地，需 AmberTools cpptraj + prod.nc）
   - 模板：上传 `an*.in` 到 `~/md/out/<system>/`，逐个跑 `cpptraj -i an*.in > an*.log 2>&1`
   - 输出回拉本地 `results/<system>/cpptraj_raw/`
2. **跑 MMPBSA**（如已有 FINAL_DECOMP_MMPBSA.dat 可跳过；建议用 gmx_MMPBSA 1.6.1 / AmberTools 24 MMPBSA.py）
3. **跑 fig01-08**：`python plot_pretty_figs.py` —— 必须先在脚本里设好 6 个体系常量（CASE / RES / OUT / D / NFR / DT_NS / T / parse_decomp）
4. **跑 fig09-22**：`python plot_pretty_figs2.py` —— 自动 import 上一脚本，复用 D 字典/风格函数
5. **跑 fig23-30**：`python plot_pretty_figs3.py` —— 同样 import 上一步；新增 8 张 DCCM/per-residue ΔG/occupancy 分布/双通道/DSSP timeline/RMSF vs B/DNA 形态/Mn RDF
6. **核图**：检查 30 张 PNG，特别关注 fig14（DSSP mapping）、fig16（DNA 残基编号偏移 = 255）、fig23（DCCM vlim）、fig28（RMSF vs B 校验）

## 图集清单（30 张）
| # | 文件 | 内容 | 数据源 |
|---|---|---|---|
| fig01 | `fig01_rmsd.png` | RMSD 时程 (蛋白+ssDNA) | `rmsd_*.dat` |
| fig02 | `fig02_rmsf.png` | RMSF (per residue) | `rmsf_prot.dat` |
| fig03 | `fig03_rg.png` | 回转半径 Rg | `rg_*.dat` |
| fig04 | `fig04_sasa.png` | SASA (cplx + 分离 buried) | `sasa_*.dat` |
| fig05 | `fig05_pca.png` | PC1-PC2 散点 + FEL | `pca_*.dat` |
| fig06 | `fig06_fel.png` | 自由能形貌图 | `pca_*.dat` |
| fig07 | `fig07_mmpbsa.png` | MM-PBSA 逐残基分解 | `FINAL_DECOMP_MMPBSA.dat` |
| fig08 | `fig08_mn_coordination.png` | 金属三齿配位距离 | `mn_*.dat` |
| fig09 | `fig09_energy_convergence.png` | 能量/温度/密度收敛 | `prod.out` |
| fig10 | `fig10_rmsd_convergence.png` | RMSD 前半 vs 后半 KDE | `rmsd_*.dat` |
| fig11 | `fig11_rmsd_rg_joint.png` | RMSD-Rg 联合分布 | `rmsd_prot.dat` + `rg_prot.dat` |
| fig12 | `fig12_pca_3d.png` | PC1-3 三维轨迹 | `pca_*.dat` |
| fig13 | `fig13_ramachandran.png` | Ramachandran 全残基 | `phi.dat` + `psi.dat` |
| fig14 | `fig14_dssp.png` | 二级结构 (DSSP) | `dssp.dat` + `dssp_sum.dat` |
| fig15 | `fig15_contacts_timeseries.png` | 界面接触数时程 | `nat_series.dat` |
| fig16 | `fig16_contact_map.png` | 残基-核苷酸接触热图 | `nat_res.dat` |
| fig17 | `fig17_hbond_occupancy.png` | 界面氢键占用率 Top 15 | `hb_all_avg.dat` |
| fig18 | `fig18_pocket_rmsd.png` | 催化口袋 RMSD | `pock_rmsd.dat` |
| fig19 | `fig19_dna_end2end.png` | DNA 首尾距 | `dna_ee.dat` |
| fig20 | `fig20_rmsf_zoom.png` | RMSF 局部放大 | `rmsf_prot.dat` |
| fig21 | `fig21_nucleotide_decomp.png` | MM-PBSA 逐核苷酸分解 | `FINAL_DECOMP_MMPBSA.dat` |
| fig22 | `fig22_hbond_timeseries.png` | 6 条关键界面氢键距离时程 | `hbd_*.dat` (an3 输出) |
| fig23 | `fig23_dccm.png` | 蛋白运动互相关矩阵 DCCM | `dccm_prot.dat` (an4b 输出) |
| fig24 | `fig24_mmpbsa_prot_residue.png` | MM-PBSA 蛋白 per-residue ΔG | `FINAL_DECOMP_MMPBSA.dat` (R 侧) |
| fig25 | `fig25_hbond_occupancy_dist.png` | 界面 H 键 occupancy 分布 + 累积 | `hb_all_avg.dat` |
| fig26 | `fig26_contacts_hbond_dual.png` | 接触 + H 键双通道时序 | `nat_series.dat` + `hb_all_series.dat` |
| fig27 | `fig27_dssp_timeline.png` | DSSP 二级结构时间演化栅栏 (1500 帧) | `dssp.dat` |
| fig28 | `fig28_rmsf_vs_bfactor.png` | MD RMSF vs 晶体 B-factor 校验 | `rmsf_prot.dat` + `bfactor_prot.dat` |
| fig29 | `fig29_dna_morphology.png` | ssDNA 形态三指标 (端距/Rg/弯曲) | `dna_ee.dat` + `rg_dna.dat` + `dna_bend.dat` |
| fig30 | `fig30_mn_rdf.png` | Mn²⁺ 配位环境径向分布 | `mn_rdf.dat` (an5 输出) |

## 风格系统（plot_pretty_figs.py 内）
- **配色常量**：`C_PROT=#2B6CB0`(深蓝), `C_DNA=#E8623C`(橙), `C_CPLX=#7E57C2`, `C_ACC=#10A37F`(翡翠), `C_PUR=#8E44AD`(紫), `C_TEAL=#3A8FB7`
- **关键残基映射 `KEY_RES`**：seq→(PDB 名, group)（F88@seq67, Y170@seq149, R253@seq232 等）
- **关键功能区**：`HELIX_CAT=(62,80)` 催化螺旋, `LOOP_MOD=(168,181)` 建模环
- **辅助函数**：`style_ax()` (去上右框+grid), `stat_box()`, `add_marginal()`, `panel_label()`, `roll_mean()`, `parse_decomp()`
- **字体**：Microsoft YaHei + DejaVu Sans（已硬编码 fallback），`U+27E8 ⟨⟩` 缺失需替换为 `<>`

## 上游 cpptraj 调用（核心踩坑）
| 坑 | 解决 |
|---|---|
| `hbond :1-254 :255-278` 报错 | 用 `,` 合并 mask → `:1-254,:255-278` |
| 合并 mask 后 series 只输出 `HBALL[UU]` 总数 | 距离法替代: an3_hbdist.in 量关键 bond 距离 |
| 5' 端 DNA 残基无 P 原子 | 用 `C5'` 起算端距 |
| isolated SASA = 0 | `strip !(:mask)` 必须同时 `parmstrip` 移除原 strip 后的冗余 parm |
| 复合物 SASA 直接 `surf :mask` 输出的是绝对值 | buried = (prot_iso − prot_in_cplx) + (dna_iso − dna_in_cplx) |
| `surf` = LCPO（快） vs `molsurf` = 数值法（慢，整批别用） | 默认 `surf` |
| cpptraj `start/stop` = 帧范围不是特征向量数 | 选模式向量用 `beg/end` |
| `projection PROJ evecs NAME :1-254@CA beg 1 end 3` 缺 mask 报错 | 必须显式 mask |
| `multidihedral` psi 从残基 2 起始 | 列索引 = seq_res − 2 |
| `parmstrip` 与 `strip` 混用冲突 | 24.5+ 不用 `parmstrip+tag` 改用纯 `strip` |

## 复用定制（每个新体系要改）
- `plot_pretty_figs.py` 顶部 6 个常量：
  ```python
  CASE = "<system_id>"      # 子目录名
  RES  = f"results/{CASE}"
  OUT  = f"{RES}/figures_pretty"
  D = dict(...)             # 把 rmsd_prot/rg_prot/rmsf_prot/sasa_* 等读入
  NFR, DT_NS = 1500, 0.1    # 1500 帧 × 0.1 ns/帧
  T = np.arange(NFR) * DT_NS
  ```
- `KEY_RES`：`{seq_idx: (pdb_name, group)}` 标注
- `HELIX_CAT / LOOP_MOD`：功能区 seq 范围
- `S1 = "<24nt 序列>"`（用于 fig16 核苷酸标注）
- 6 个 cpptraj 脚本里 parm/trajin 路径按需改

## 已验证坑清单（写新体系前必读）
- **DSSP code 硬编码**：`{0:coil, 1:E, 2:B, 3:G, 4:H, 5:I, 6:T, 7:S}`；`dssp_sum.dat` 第 0 列是 `#Residue` 编号，**必须先 `[:,1:]`**
- **nat_res.dat DNA 编号偏移**：parm 255=nt1（不是 nt0），偏移常量 254，不要用 `dna_min-1`
- **hb_all_avg.dat 过滤界面**：`is_if = (ra<=254<rd) or (rd<=254<ra)`，因 seq 蛋白 1-254 而 DNA 255-278
- **Microsoft YaHei 缺 `⟨⟩ U+27E8/27E9`**：用 `<>` 或 `mean(...)` 文本替代
- **`pLDDT`/`ipTM` 同平台比较**：Boltz vs AF3 系统偏差 0.05-0.1
- **Kabsch 铁律**：跨预测比 RMSD 必须先叠加（1800 骨架=1 个的教训）
- **Boltz ModelCIF** 用 `parts[3]` 原子名（`parts[2]` 元素符号 bug 曾致距离高估）

## 已知未解决
- ~~**fig22 完整时序版**（6 条关键氢键距离演化）~~ ✅ **已解决（2026-09-02）**：
  - CHPC 跑 `an3_hbdist.in` → 6 个 `hbd_*.dat`（1500 帧 × 0.1 ns/帧）
  - fig22 已重写为「上=6 键逐帧距离时序（细线+5 ns 滚动均值+3.5 Å 阈值）；下=成键状态栅栏 (binary raster)」
  - 跑法见下方"CHPC cpptraj 加载"速查

## CHPC cpptraj 加载速查（最常踩坑）
CHPC 登录节点 `workstation`，`which cpptraj` 为空。amber/24 自带 cpptraj 链接了系统不提供的 `libopenblas/libplumed`，spack ubuntu 版运行库与系统 centos7 glibc 不兼容。**唯一可用的现场方案**：
```bash
G10=/opt/app/spack/opt/spack/linux-centos7-haswell/gcc-4.8.5/gcc-10.1.0-2new4oxsi6o5ejrrjxsjtvvxubujfyyk
export LD_LIBRARY_PATH=$G10/lib64:/opt/app/amber22_gpu/lib
/opt/app/amber22_gpu/bin/cpptraj --version     # → CPPTRAJ V6.4.4
```
- 绝不要用 `/opt/app/amber/24/bin/cpptraj`（报 libarpack/libopenblas/libplumed/libgfortran 缺）
- 绝不要用 `/opt/app/amber22_gpu/bin/cpptraj` 不带 LD_LIBRARY_PATH（报 libgfortran.so.5 缺）
- prod.nc 1.3 GB **不要拉回本地**；cpptraj 算后只拉 .dat
- 登录节点跑轻量 cpptraj OK（1500 帧 distance < 30 s），仅 GPU / 长作业需 sbatch
- 凭据：`CHPC_PASS=...`（paramiko 5）；推荐 `md_specificity/ssh_run.py` 一键执行/拉取/推送
- ensemble-scorer：未经 PprI 实测校准，排序可信绝对值无精度承诺

## 调用方式（其他项目）
1. 把 `scripts/plot_pretty_figs.py` + `plot_pretty_figs2.py` + `plot_pretty_figs3.py` 拷到新项目
2. 把 `cpptraj/an*.in` 拷到新项目，按需修改 parm/trajin 路径
3. 改 6 个常量（CASE/RES/OUT/D/NFR/DT_NS/T）和 KEY_RES/S1
4. 跑 cpptraj → 拉数据 → `python plot_pretty_figs.py` → `python plot_pretty_figs2.py` → `python plot_pretty_figs3.py`
5. 跑 HTML 图廊（30 张）：`python scripts/figures_gallery.py` → `examples/index.html`

## GitHub 仓库
- 公开仓库 [`Elephenman/md-amber-figure-pipeline`](https://github.com/Elephenman/md-amber-figure-pipeline)
- 收录本 skill 全部文件 + 30 张示例 PNG + 深色科研风 HTML 图廊 + 8 个 cpptraj 输入
- 接受 PR/issue：欢迎补充你自己的新体系模板
