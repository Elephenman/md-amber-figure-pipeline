# 流水线架构与集群部署 — PIPELINE_ARCHITECTURE

> 本文件说明 `md_easy.sh` 的阶段架构、产物布局、以及如何在集群（CHPC 登录节点 / sbatch GPU 分区）上跑完整条流水线。
> 核心设计目标：**一条命令、断点续跑、傻瓜式交互**；仓库可整体拷贝到集群，不改任何路径即可跑。

---

## 1. 总体架构

```
输入 PDB (蛋白+核酸[+Mn])
      │  bash md_easy.sh   （唯一入口）
      ▼
┌──────────────────────────────────────────────────────────────┐
│ md_easy.sh —— 参数解析(交互/配置) + 断点调度(0-6 阶段循环)      │
└──────────────────────────────────────────────────────────────┘
      │ source 公共库
      ▼
pipeline/common.sh —— 仓库定位/配色日志/默认配置/派生路径/引擎探测/命令执行(dry-run 感知)
      │
      ▼  逐阶段执行（每阶段成功写 .md_easy/<SYSTEM>/stageN.done）
┌──────┬──────────┬──────────────┬──────────┬──────────┬─────────┐
│stage0│ stage1   │ stage2       │ stage3   │ stage4   │ stage5  │ stage6
│ 体检 │ tleap    │ 约束/配位     │ MD       │ cpptraj  │ MM-PBSA │ 出图
└──────┴──────────┴──────────────┴──────────┴──────────┴─────────┘
```

| 层 | 内容 |
|---|---|
| **总入口** | `md_easy.sh`（唯一需要敲的命令；交互问答 / `-c config` 非交互） |
| **公共库** | `pipeline/common.sh`：加载默认值、路径派生、引擎探测、`render()` 模板渲染、`exec_cmd()` dry-run 感知执行 |
| **阶段脚本** | `pipeline/stage0..6_*.sh`：每个阶段独立可跑（自动 source common） |
| **模板** | `templates/*.tmpl`：leap.in / min1 / min2 / heat / equil / prod / mmpbsa.in |
| **cpptraj 模板** | `cpptraj/an*.in`：an0 基础量在 stage4 内联拼装；an1..an5 用模板渲染 |
| **工具** | `tools/*.py`：check_input_pdb / fix_his / make_restraints / render_tmpl / ssh_run |
| **绘图** | `scripts/plot_pretty_figs{1,2,3}.py` + `figures_gallery.py` |
| **力场** | `ff/mn_cm12-6.frcmod`（12-6 Mn²⁺）+ `.itp`（GROMACS 备用） |

## 2. 配置与断点机制

- **配置三层**：`config/project.env.example`（文档注释最全）→ 交互问答生成的 `project.env` → 运行中的环境变量默认值兜底（`load_defaults`）。
- 最小必需项：`INPUT_PDB` + `SYSTEM`；其余全有默认。
- 断点目录：`.md_easy/<SYSTEM>/stage{N}_{name}.done`；重跑跳过已完成阶段。
  - `-f` 强制全跑；`--from N` 从 N 跑起（自动清 ≥N 标记）；`-s 1,3-5` 选跑；`--dry-run` 只体检预览；`--env-only` 只校验配置。
- **单阶段脚本可独立执行**：`bash pipeline/stage6_figs.sh`（需要环境变量，建议仍走 md_easy.sh）。

## 3. 产物布局

```
<RUN_DIR>/                 # 默认 = 调用 md_easy.sh 的目录
├── project.env            # 交互生成的配置
├── .md_easy/<SYSTEM>/     # 断点标记
├── logs/<SYSTEM>/         # 阶段日志
├── topol/<SYSTEM>/        # v.parm7, v.rst7, input_hid.pdb, masks.env, disang.txt, leap.log
├── out/<SYSTEM>/rep{N}/   # min1/2, heat, equil, prod 的 {in,out,rst,nc,mdinfo}
└── results/<SYSTEM>/      # 分析数据：*.dat(顶层) + cpptraj_raw/ + prod.out + figures_pretty/
```

## 4. 引擎 / 平台策略

- `ENGINE=auto`：`detect_amber()` 探测 `AMBERHOME/bin/pmemd.cuda`（GPU）→ 有则用；任何阶段 GPU 失败自动退 `pmemd`（CPU）重跑该步（stage3 内 `|| { warn; ...PMEMD_CPU...; }`）。
- stage0 额外校验 python + numpy/scipy/matplotlib（出图用）与 nvidia-smi。
- **Windows (Git Bash)**：仓库已做 MSYS 路径归一化（脚本内 `/a/...` → `A:/...`），`.gitattributes` 强制 LF 行尾，可直接在本地 Git Bash 跑 stage4/5/6 的 CPU 部分。

## 5. CHPC 集群部署（浙大 10.202.94.52:20009，实测）

### 5.1 登录节点环境

```bash
# cpptraj 无 PATH —— 显式加载（V6.4.4）：
export LD_LIBRARY_PATH=/opt/app/spack/opt/spack/linux-centos7-haswell/gcc-4.8.5/gcc-10.1.0-2new4oxsi6o5ejrrjxsjtvvxubujfyyk/lib64:/opt/app/amber22_gpu/lib
#   ↑ spack centos7 gcc-10 lib64 提供 libgfortran.so.5 + 新 libstdc++
CPPTRAJ=/opt/app/amber22_gpu/bin/cpptraj     # 实测可用
# ⚠️ 绝不要 /opt/app/amber/24/bin/cpptraj（缺 libopenblas/libplumed；spack ubuntu 库与系统 centos7 glibc 不兼容）
```

### 5.2 一条命令跑通建议路径

登录节点只跑**轻量 CPU 任务**（cpptraj 距离/接触 <30 s 量级没问题）。全流程推荐分三段：

```bash
# ① 拷贝仓库 + 输入 PDB 到集群
git clone https://github.com/Elephenman/md-amber-figure-pipeline.git   # 或 scp 整个目录
cd md-amber-figure-pipeline && cp /path/to/complex.pdb .
bash md_easy.sh                 # 交互填 PDB → 生成 project.env（或手工复制 config/project.env.example）

# ② 登录节点建系（0-2 阶段，CPU 分钟级；3 阶段提交 GPU）
bash md_easy.sh -c project.env -s 0,1,2
# ③ GPU sbatch 跑 prod（见 5.3）→ 回登录节点跑后处理
bash md_easy.sh -c project.env --from 4    # cpptraj 采集 + MM-PBSA + 30 图
```

### 5.3 sbatch GPU 提交模板（4090 分区实测参数）

> 分区 `4090`（4090D）/ `gpu`（V100）；**禁 `--mem`**；计费注释 `--comment=ls_lhz`；提交前 `mkdir -p logs`。

```bash
#!/bin/bash
#SBATCH --job-name=md_prod
#SBATCH --partition=4090
#SBATCH --output=logs/md_prod_%j.out
#SBATCH --error=logs/md_prod_%j.err
#SBATCH --comment=ls_lhz
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
# 加载 amber（GPU 版）
export AMBERHOME=/opt/app/amber22_gpu
export PATH=$AMBERHOME/bin:$PATH
export LD_LIBRARY_PATH=$AMBERHOME/lib:$LD_LIBRARY_PATH
cd $SLURM_SUBMIT_DIR
# 只跑阶段 3（min/heat/equil/prod 全在 GPU），断点续跑由脚本自身保证
bash md_easy.sh -c project.env -s 3
```

**集群踩坑备忘（项目实战验证）**：
- sbatch 脚本**严禁放入 md_easy 的输入/输出目录**，放仓库根目录同级；`.sbatch` 放错位置会触发 parse error。
- 提交前 `mkdir -p logs`，否则 sbatch 报错找不到输出目录。
- sbatch 内含 `$` 变量（如 `$SLURM_SUBMIT_DIR`）时**不要**用 `--export` 内联传参；直接写在脚本内最稳。
- 上传文本文件必须 LF 行尾（Windows CRLF 会让 bash/cpptraj 报 `\r` 错误）—— 用 `sed -i 's/\r$//'` 或 git clone。
- **大轨迹（prod.nc 常 >1 GB）勿拉回本地**：登录节点跑 cpptraj 只拉 `.dat`。

### 5.4 登录节点一键后处理（轻量版）

```bash
# 只做数据采集 + 出图（CPU，cpptraj 单帧量级任务）
bash md_easy.sh -c project.env --from 4
# 若 prod.nc 在集群、想本地出图：只拉 results/ 与 topol/<SYSTEM>/{v.parm7,masks.env}，再本地跑
bash md_easy.sh -c project.env -s 6        # 本地 stage6（脚本内置 MSYS 路径归一化）
```

## 6. 远程工具 tools/ssh_run.py（paramiko）

```bash
# 依赖: python3 + paramiko（pip install paramiko）
CHPC_PASS=xxx python tools/ssh_run.py "ls ~/md_amber_figure_pipeline"   # 执行远程命令
CHPC_PASS=xxx python tools/ssh_run.py pull /path/prod.nc                 # 拉文件到 ./mirror/
CHPC_PASS=xxx python tools/ssh_run.py push local.pdb /remote/abs/path    # 推送本地文件
```

## 7. 阶段失败速查

| 现象 | 根因 | 处理 |
|---|---|---|
| stage0 die: 缺 AmberTools | AMBERHOME 未设/不完整 | `source amber.sh` 或 `export AMBERHOME=...` |
| stage1 die: Mn 电荷≠+2 | frcmod 未加载/PDB Mn 电荷列错误 | 查 leap.log；核对 `ff/mn_cm12-6.frcmod` |
| stage3 prod RESTRAINT>100 | 起始配位几何崩坏 | 见 `docs/METAL_CHECKLIST.md` §2/§5 |
| stage4 die: prod.nc 不存在 | 阶段 3 未跑完/被 kill | `squeue` 查作业；重提交 `-s 3`（自动续跑） |
| stage6 die: 缺数据文件 | 数据未齐 | 提示会列出缺哪个文件；`--from 4` 补跑 |
| 绘图缺可选数据（fig07/08/21/22/24/28/30） | MM-PBSA 未跑 / 无金属 / 无 bfactor | 脚本自动占位图，不中断 |
