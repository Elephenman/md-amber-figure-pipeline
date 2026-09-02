---
name: md-amber-figure-pipeline
description: 端到端 Amber MD 全流程自动化流水线——从输入 PDB 一条命令跑到 30 张推文级图 + HTML 画廊（环境体检→tleap 建拓扑→Mn2+/DISANG 约束→min/heat/equil/prod→cpptraj 采集→MM-PBSA→30 图）。傻瓜式交互，无需调参。适用于蛋白-核酸(-金属)复合物 MD 从头建系到论文图全套。触发词：MD 漂亮图、MD 美化作图、Amber 流水线、建拓扑跑 MD、一条命令跑 MD、md_easy、cpptraj 出图、30 图集、WT MD 出图、150 ns MD。
agent_created: true
---

# md-amber-figure-pipeline — Amber MD 一键全流程（输入 PDB → 30 图）

## 一句话
把一份「蛋白(+核酸)(+金属)复合物 PDB」喂给 `md_easy.sh`，自动完成 Amber MD 全流程并产出 **30 张推文级 PNG + HTML 画廊**。全部默认参数来自 PprI·ssDNA·Mn²⁺ 150 ns 真实项目验证，**傻瓜式无需调参**。

## 快用（3 条命令）
```bash
bash md_easy.sh                        # 交互问答（回车=默认），生成 project.env
bash md_easy.sh -c my.env              # 正式跑全流程（0体检→1tleap→2约束→3MD→4cpptraj→5MMPBSA→6图）
bash md_easy.sh -c my.env --from 4     # 断点续跑（从阶段 4）
# 产物: topol|out|results/$SYSTEM/figures_pretty/fig01..30.png + index.html
```

## 适用对象与前提
- 蛋白·核酸·金属复合物（PprI 模板：蛋白 254aa + 24nt ssDNA + Mn²⁺，~30-100k 原子，1500 帧×0.1ns）
- **前提**：装了 AmberTools（tleap/pmemd/cpptraj/MMPBSA.py）；有 NVIDIA GPU 可加速（CPU 自动兜底）
- 残基编号：蛋白 seq 1..N，DNA 从 N+1 起（PprI: parm 255=nt1）；跨体系由 make_restraints 自动输出 masks.env

## 阶段架构
| 阶段 | 脚本 | 干什么 |
|---|---|---|
| 0 | stage0_env.sh | 引擎探测 + PDB 体检 + **配位 His→HID/HIE 自动改写**（输出 topol/input_hid.pdb，不污染原输入） |
| 1 | stage1_tleap.sh | tleap: ff19SB/OL15/TIP3P + Na+/Cl- 中和 + 金属 frcmod |
| 2 | stage2_restraints.sh | Mn 配位 DISANG(masks.env/disang.txt) + 配位距离体检 |
| 3 | stage3_md.sh | min1(wt25)→min2→heat(NVT wt10)→equil(NPT wt2)→**prod(NPT ntr=0 自由)**，GPU 失败自动退 CPU |
| 4 | stage4_cpptraj.sh | an0 基础量 + an1..an5 专项（体系类型自适应跳过无关分析） |
| 5 | stage5_mmpbsa.sh | MM-PBSA 逐残基分解（可选 DO_MMPBSA=yes，需 MMPBSA.py） |
| 6 | stage6_figs.sh | 30 图 + 画廊（数据预检硬失败、可选数据占位图） |

## 30 图数据源速查（fig 编号 ↔ 文件）
- fig01-06（RMSD/RMSF/Rg/SASA/PCA/FEL）← results 顶层 rmsd/rmsf/rg/sasa/pca_*.dat
- fig07（MM-PBSA 能量项）← FINAL_RESULTS.dat+FINAL_DECOMP_MMPBSA.dat（缺→占位）
- fig08（Mn 三齿配位）← mn_H71/H75/E102a/E102b.dat（无金属→占位）
- fig09（能量收敛）← prod.out；fig10-12 收敛/联合/3D ← 顶层；fig13-14 ← cpptraj_raw phi/psi/dssp
- fig15-17（接触/hbond）← an2b nat_*/hb_all_*；fig18-19 ← pock_rmsd/dna_ee；fig20 ← rmsf
- fig21（逐核苷酸）← FINAL_DECOMP（缺→占位）；fig22（6 键距离）← hbd_*.dat（config HBOND_PAIRS，缺→占位）
- fig23（DCCM）← an4b dccm_prot.dat；fig24（蛋白 per-residue）← FINAL_DECOMP（缺→占位）
- fig25-26 ← an2b；fig27 ← dssp.dat；fig28 ← bfactor_prot.dat（**可选输入**：`#Seq PDBres Bfactor_CA`，缺→占位）
- fig29 ← dna_ee/rg_dna/dna_bend；fig30 ← mn_rdf.dat（无金属→占位）

## 关键默认参数（=项目验证结论，勿轻易改）
| 参数 | 默认 | 铁律/依据 |
|---|---|---|
| PROT_FF/DNA_FF/WAT | ff19SB/OL15/TIP3P | Amber 蛋白·核酸标准 |
| 约束 | heat wt10 → equil wt2 → **prod ntr=0** | **prod 骨架约束(哪怕 wt1.0)会把 RMSD 钉死 0.39Å → 病态采样，必须自由** |
| Mn²⁺ | 12-6 frcmod(ff/mn_cm12-6.frcmod) + flat-bottom DISANG(2.1-2.5Å rk=10) | sander GB 不识别 Mn 类型 → MMPBSA strip_mask 剥 Mn |
| 配位 His | Mn 3.0Å 内 HIS→HID(ε去质子)/HIE(δ去质子) | 金属配位 His 必须对质子化态正确 |
| HBOND_PAIRS | config 定义「名字=残基:原子=残基:原子」`\|` 分隔 | fig22 数据源，PprI 模板 6 键见 config 注释 |

## 体系类型支持
| 体系 | stage4 | stage6 |
|---|---|---|
| 蛋白+DNA+金属 | 全量 | 30 图 |
| 蛋白+DNA 无金属 | 跳过 Mn 分析 | 30 图（fig08/30 占位） |
| apo 裸蛋白 | 跳过 DNA/Mn | 跳过出图（图集面向复合物） |

plot1 (fig01-08) 已通用化：`MDEASY_RES` env 覆盖数据目录 + MSYS 路径归一化（Git Bash 下 /a/... → A:/...）+ `_ttl()` 标题前缀用 CASE。plot2/3 (fig09-30) 为 PprI 案例深度绑定（注释/suptitle），复用新体系时按需裁剪或仅用 plot1。

## 断点 / 排查
- 断点文件 `.md_easy/$SYSTEM/stageN.done`；`-f` 强制、`--from N` 续跑、`-s` 单跑
- `--dry-run` 只体检+预览；`--env-only` 只校验配置
- 金属体系排查顺序：① stage0 是否改写 HID/HIE → ② input_hid.pdb 里 Mn-N/O 1.9-2.6Å → ③ frcmod 存在 → ④ masks.env 的 COORD1..4
- 常见失败：stage6 缺数据 → 提示缺哪个文件 + `--from 4` 补跑；绘图缺可选数据 → 自动占位图（fig07/08/21/22/24/28/30）

## 集群速查（CHPC 10.202.94.52:20009 登录节点）
- cpptraj 无 PATH：需 `export LD_LIBRARY_PATH=/opt/app/spack/.../gcc-10.1.0-.../lib64:/opt/app/amber22_gpu/lib` + `/opt/app/amber22_gpu/bin/cpptraj`（V6.4.4）
- **绝不要** `/opt/app/amber/24/bin/cpptraj`（缺 libopenblas/libplumed；spack ubuntu 库与 centos7 glibc 不兼容）
- 提交 sbatch 前 `mkdir -p logs`；sbatch 禁入 boltz 输入目录等坑见 docs/PIPELINE_ARCHITECTURE.md
- 大轨迹（>1GB）勿拉回本地，cpptraj 算后只拉 .dat；远程工具 tools/ssh_run.py（paramiko）

## 已知坑（继承自项目实战）
- DSSP code 硬编码；`dssp_sum.dat` 先 `[:,1:]`
- nat_res DNA 偏移 = 蛋白残基数（parm 255=nt1）
- hbond 双 mask 用 `,` 合并 `:1-254,:255-278`；合并后 series 只有 HBALL[UU] 总数 → 逐键用 an3_hbdist 距离法
- 孤立 SASA 必须 strip+parmstrip 同步；DCCM 必须先 rms fit 去刚体
- Microsoft YaHei 缺 `⟨⟩` 字形 → 用 `<>`
- Windows Git Bash 跑 python：脚本已内置 MSYS 路径归一化
- 出图需要 numpy/scipy/matplotlib（无 GPU 需求）

## 调用方式（其他项目复用）
1. 复制整个仓库或子集（md_easy.sh + pipeline/ + templates/ + cpptraj/ + tools/ + scripts/ + ff/ + config/）
2. 跑 `bash md_easy.sh` 交互填 PDB → 或复制 `config/project.env.example` 改 INPUT_PDB/SYSTEM
3. 若要 fig22：在 project.env 填 HBOND_PAIRS（PprI 模板 6 键格式见 config 注释）
4. 若要 fig28：放 bfactor_prot.dat 到 results/<SYSTEM>/（可选）
5. 收图：results/<SYSTEM>/figures_pretty/ + index.html

## GitHub 仓库
- 公开仓库 [`Elephenman/md-amber-figure-pipeline`](https://github.com/Elephenman/md-amber-figure-pipeline)
- 收录全部流水线脚本 + 模板 + 30 张示例 PNG + 深色科研风画廊 + docs（PIPELINE_ARCHITECTURE / METAL_CHECKLIST / METHODS_TEMPLATE）
- 接受 PR/issue：欢迎补充新体系模板
