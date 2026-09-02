# md-amber-figure-pipeline

**端到端 Amber MD 自动化流水线**：一条命令完成 输入 PDB → 环境体检 → tleap 建拓扑 → Mn²⁺/约束方案 → min/heat/equil/prod → cpptraj 数据采集 → MM-PBSA(可选) → 30 张推文级图 + HTML 画廊。

> 模板与默认参数全部来自 **PprI·ssDNA·Mn²⁺ 150 ns MD** 真实项目验证（[Elephenman/PprI_ssDNA_design](https://github.com/Elephenman/PprI_ssDNA_design)）。仓库只装脚本/模板/示例输出，不装轨迹数据。

---

## 三步傻瓜使用

```bash
# ① 一键交互（全部回车=默认，自动探测金属/体系类型，生成 project.env）
bash md_easy.sh

# ② 正式模拟：配置跑全流程（0 体检 → 1 tleap → 2 约束 → 3 MD → 4 cpptraj → 5 MM-PBSA → 6 图）
bash md_easy.sh -c my.env

# ③ 看结果
#    30 张图: results/<SYSTEM>/figures_pretty/  （打开 index.html 画廊）
#    断点续跑: bash md_easy.sh -c my.env --from 3   （从阶段 3 继续）
```

**不需要调参**：默认参数即项目验证结论（见下），只填输入 PDB 即可。

| 命令 | 作用 |
|---|---|
| `bash md_easy.sh` | 交互问答生成 `project.env`（回车直通） |
| `bash md_easy.sh -c my.env` | 非交互跑全流程（正式模拟推荐） |
| `bash md_easy.sh -c my.env -f` | 强制重跑（忽略断点） |
| `bash md_easy.sh -c my.env --from 3` | 从阶段 3 续跑到结束 |
| `bash md_easy.sh -c my.env -s 1,3-5` | 只跑指定阶段（逗号/区间） |
| `bash md_easy.sh -c my.env --dry-run` | 只体检 + 预览命令，不执行 |
| `bash md_easy.sh -c my.env --env-only` | 只校验配置后退出 |

## 目录分类

```
md-amber-figure-pipeline/
├── md_easy.sh                 # ★ 傻瓜式总入口（唯一需要敲的命令）
├── config/project.env.example # 全部可调参数(全注解) —— 复制为 project.env 或交互生成
├── pipeline/                  # 6 个阶段脚本 + 公共库
│   ├── common.sh              #   路径/日志/渲染/断点/引擎探测 公共函数
│   ├── stage0_env.sh          #   环境体检 + PDB 体检 + His 质子化态改写(HID/HIE)
│   ├── stage1_tleap.sh        #   tleap 拓扑 (ff19SB/OL15/TIP3P + Na+/Cl- + 金属 frcmod)
│   ├── stage2_restraints.sh   #   Mn 配位 DISANG + 掩码 (masks.env)
│   ├── stage3_md.sh           #   min1→min2→heat→equil→prod (GPU 自动退 CPU)
│   ├── stage4_cpptraj.sh      #   an0 基础量 + an1..an5 专项量采集
│   ├── stage5_mmpbsa.sh       #   MM-PBSA 逐残基分解 (可选, DO_MMPBSA=yes)
│   └── stage6_figs.sh         #   30 图一键 + HTML 画廊
├── templates/                 # leap.in / min/heat/equil/prod.in / mmpbsa.in 模板
├── cpptraj/                   # an1_struct/an2b_contacts/an2c_pock/an2d_ee/
│                              #   an3_hbdist/an4_dynamics/an4b_dccm_fit/an5_rdf 模板
├── tools/                     # check_input_pdb / fix_his / make_restraints /
│                              #   render_tmpl / ssh_run
├── scripts/                   # plot_pretty_figs.py(01-08) / 2(09-22) / 3(23-30) /
│                              #   figures_gallery.py(通用画廊)
├── ff/                        # mn_cm12-6.frcmod (Mn²⁺ 12-6 非键参数)
└── examples/                  # 30 张示例 PNG（PprI·WT·S1 150 ns 真实 MD）
```

## 产物布局（自动创建）

```
$PWD/topol/<SYSTEM>/        # v.parm7/v.rst7/input_hid.pdb/masks.env/disang.txt
$PWD/out/<SYSTEM>/rep{N}/   # min/heat/equil/prod 轨迹 + 日志
$PWD/results/<SYSTEM>/      # 分析数据: *.dat(顶层) + cpptraj_raw/ + prod.out
$PWD/results/<SYSTEM>/figures_pretty/   # fig01..fig30.png + index.html
```

## 30 张图总览

| 图号 | 内容 | 数据源 |
|---|---|---|
| fig01-04 | RMSD / RMSF / Rg / SASA(cplx vs isolated → buried) | an0 顶层 .dat |
| fig05-06 | PCA + FEL 自由能形貌 | an0 pca |
| fig07 | MM-PBSA GB/PB 能量项 + ΔGbind | stage5 (缺→占位) |
| fig08 | Mn²⁺ 三齿配位距离时序 (H92/H96/E123) | an0 mn_*.dat (无金属→占位) |
| fig09-12 | 能量收敛 / RMSD 前后半 KDE / RMSD-Rg 联合 / PC1-3 3D | prod.out + 顶层 |
| fig13-14 | Ramachandran / DSSP 二级结构组成 | an1 |
| fig15-17 | 界面接触数 / 接触热图 / H 键 occupancy Top15 | an2b |
| fig18-20 | Mn 催化口袋 RMSD / ssDNA 首尾距 / RMSF 局部放大 | an2c/an2d |
| fig21 | MM-PBSA 逐核苷酸分解 | stage5 (缺→占位) |
| fig22 | 界面关键 H 键逐帧距离 + 成键栅栏 (config HBOND_PAIRS) | an3 (缺→占位) |
| fig23 | DCCM 蛋白运动互相关矩阵 | an4b |
| fig24 | MM-PBSA 蛋白 per-residue ΔG | stage5 (缺→占位) |
| fig25-26 | H 键 occupancy 分布 / 接触+H 键双通道 | an2b |
| fig27 | DSSP 时间演化栅栏 | an1 |
| fig28 | MD RMSF vs 晶体 B-factor 校验 | bfactor_prot.dat (可选→占位) |
| fig29 | ssDNA 形态三指标: 端距/Rg/弯曲角 | an2d/an4 |
| fig30 | Mn²⁺ 配位 RDF | an5 (无金属→占位) |

## 关键默认参数（= 项目验证结论，勿轻易改）

| 参数 | 默认 | 依据 |
|---|---|---|
| 蛋白力场 | ff19SB (可 ff14SB) | Amber 蛋白默认 |
| 核酸力场 | OL15 (可 bsc1) | ssDNA 常用 |
| 水/离子 | TIP3P + Na+/Cl- 12-6-4 中和 | Amber 默认 |
| 约束方案 | min wt=25 → heat NVT wt=10 → equil NPT wt=2 → **prod ntr=0 自由** | 见下「prod 铁律」 |
| Mn²⁺ | 12-6 非键 (ff/mn_cm12-6.frcmod) + **flat-bottom DISANG(2.1-2.5 Å, rk=10)** 锚定配位 | 见 docs/METAL_CHECKLIST.md |
| 配位 His | 自动探测 Mn 3.0 Å 内 HIS → 改写 HID(ε 去质子)/HIE(δ 去质子) | stage0 |
| 生产 | 300 K NPT, 每重复 PROD_NS ns (默认 150), 每 100 ps 一帧 | — |

### ⚠️ prod 铁律（2026-09-01 项目教训）

**prod 阶段必须 ntr=0 自由动力学**。骨架约束(即使 wt=1.0)会把蛋白 RMSD 钉死在 ~0.39 Å，产生病态采样；只保留 Mn 配位 flat-bottom DISANG 即可（模板已固化）。

### 体系类型支持

| 体系 | stage0-3 | stage4 | stage5 | stage6 |
|---|---|---|---|---|
| 蛋白+DNA+金属 (PprI 类) | ✅ | ✅ 全量 | ✅ | ✅ 30 图 |
| 蛋白+DNA (无金属) | ✅ | ✅ (跳过 Mn 分析) | ✅ | ✅ (fig08/30 占位) |
| apo 裸蛋白 | ✅ | ✅ (跳过 DNA/Mn 分析) | ⏭ 需核酸 | ⏭ 图集面向复合物 |
| freedna (纯核酸) | ✅ | ✅ (跳过蛋白分析) | 需蛋白 | ⏭ 同上 |

plot2/plot3 (fig09-30) 为 **PprI 案例深度绑定**脚本（残基/序列注释、suptitle 保留案例名）；plot1 (fig01-08) 已通用化（`MDEASY_RES` 覆盖数据目录 + `_ttl()` 标题前缀用 CASE）。新蛋白体系跑图：推荐复用 config 与 plot1，或用 `DO_FIGS=no` 跳过出图。

## 常用命令速查

```bash
# 只跑后处理（假设已有 topol/out）
bash md_easy.sh -c my.env --from 4          # cpptraj → MM-PBSA → 图
bash md_easy.sh -c my.env -s 6              # 只出图（数据已齐）

# 金属体系检查清单（配位 His/Glu 探测失败时的排查顺序）
#   1) stage0 fix_his 是否把配位 His 改写为 HID/HIE
#   2) topol/<SYSTEM>/input_hid.pdb 里 Mn 与 N/O 距离 1.9-2.6 Å
#   3) FRCMOD 是否存在（ff/mn_cm12-6.frcmod）
# 详见 docs/METAL_CHECKLIST.md
```

## 已知坑（必读）

- **Windows (Git Bash) 跑 python 出图**：脚本已内置 MSYS 路径归一化（`/a/...` → `A:/...`），无需手工转换。
- **prod 必须 ntr=0**（见上）。
- **DSSP code 硬编码** `{0:coil,1:E,2:B,3:G,4:H,5:I,6:T,7:S}`；`dssp_sum.dat` 先 `[:,1:]`。
- **nat_res.dat DNA 偏移**：默认蛋白残基 + 1 即 nt1（PprI: parm 255=nt1）；跨体系需核对 make_restraints 输出的 DNA_RANGE。
- **hbond 双 mask 必须 `,` 合并** `:1-254,:255-278`；合并后 series 只有 `HBALL[UU]` 总数 → 逐键距离用 an3_hbdist.in。
- **孤立 SASA**：`strip + parmstrip` 同步剥除，否则结果恒 0。
- **DCCM** 前必须 `rms fit first` 去刚体，否则全正相关（刚体扩散主导）。
- **MMPBSA strip_mask 剥 Mn**：sander GB 不识别 Mn 原子类型 → rec/lig/cp 拓扑与轨迹都要剥 Mn（stage5 已自动处理）。
- **Microsoft YaHei 缺 `⟨⟩` 字形**：图中用 `<>` 或文字替代。

## 集群 / 引擎

- `ENGINE=auto`：自动探测 `pmemd.cuda`(GPU) → 失败自动退 `pmemd`(CPU)。
- 集群 sbatch 提交模板与登录节点一键依赖链见 `docs/PIPELINE_ARCHITECTURE.md`。
- 远程收数据：`tools/ssh_run.py`（paramiko 封装，支持 cmd/pull/push）。

## 依赖

- AmberTools 23+（tleap / pmemd / cpptraj / MMPBSA.py）—— 跑 MD 必需
- Python 3.10+：`numpy scipy matplotlib`（出图用，无 GPU 需求）
- 可选 NVIDIA GPU（pmemd.cuda 加速，CPU 可兜底）

## 复现示例

`examples/` 下 30 张 PNG 是 **PprI(WT)·S1-ssDNA·Mn²⁺** 150 ns 真实 MD 输出（蛋白 seq 1-254，ssDNA nt 1-24 = `TCATGAGCAGTTTTTTGTTTTTTT`）。数值摘要：
蛋白 RMSD 1.93±0.36 Å · Mn 三齿配位 H92/H96 严格保持 · 界面 H 键 R253–T24 76.7% · MM-PBSA ΣΔG −126.2 kcal/mol · R207 单残基 −13.0 top1 · MD RMSF vs 晶体 B Pearson r=0.62。

## 论文 Methods 模板

含 AMBER 选型辩护（为什么 AMBER 而非 GROMACS）、金属 12-6 参数化、约束协议英文段落 → `docs/METHODS_TEMPLATE.md`。

## License

MIT
