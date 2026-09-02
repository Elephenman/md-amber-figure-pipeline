# 金属体系检查清单 — Mn²⁺ 配位型复合物 MD

> 适用范围：PDB 含 `MN` HETATM 的蛋白·核酸·金属复合物（PprI 模板：HEXXH 基序 H92/H96/E123 三齿配位 Mn²⁺）。
> 流水线自动处理大部分问题，本清单用于**失败时按序排查**与**运行前人工预检**。
> 相关代码：`stage0_env.sh`（His 改写）→ `stage1_tleap.sh`（frcmod）→ `stage2_restraints.sh`（DISANG+配位体检）→ `stage3_md.sh`（RESTRAINT 能量闸门）。

---

## 0. 快速排查顺序（代码内提示的入口）

| 步骤 | 查什么 | 失败信号 | 见 |
|---|---|---|---|
| ① | stage0 是否把配位 His 改写为 HID/HIE | `METAL=no` 或 fix_his 未执行 | §1 |
| ② | `topol/<SYSTEM>/input_hid.pdb` 中 Mn 与配位 N/O 距离 | 距离 >3.5 Å | §2 |
| ③ | frcmod 是否存在并被加载 | `die: 金属 frcmod 不存在` | §3 |
| ④ | `masks.env` 的 COORD1..4 / MN_RES 是否合理 | MN_RES=- 或 COORD 缺失 | §4 |
| ⑤ | prod 首帧 RESTRAINT 能量 | >100 kcal/mol（stage3 直接报错退出） | §5 |

---

## 1. 金属判定与 His 质子化态

- `METAL=auto`：stage0 用 `grep '^HETATM.* MN '` 探测 PDB 是否含 Mn。**探测失败最常见原因：PDB 中 Mn 元素列写法不合规**（`MN` 前无两个空格等）。手动确认：`grep -n "^HETATM" 输入.pdb | grep " MN "`。
- 探测成功 → `METAL=yes`，stage0 调用 `tools/fix_his.py`：
  - 探测 Mn 3.0 Å 内的 HIS 残基，按最近配位原子改写：**ε-去质子 = HID（配位 NE2 时），δ-去质子 = HIE（配位 ND1 时）**。
  - 输出 `topol/<SYSTEM>/input_hid.pdb`（**不污染用户输入 PDB**）；stage1 自动优先用它建拓扑。
- 若配位残基**不是 His**（如 Glu/Asp/水），fix_his 不改写，但 stage2 的 make_restraints 仍会探测 Glu 羧基 O 配位并生成 DISANG —— 此时需核对 §4。

## 2. 起始结构配位几何

- stage2 用 cpptraj 对 `v.rst7` 量 4 条配位距离，打印 `c1..c4`；**健康区间 1.9–2.6 Å**。
- 若 >3.5 Å：说明输入 PDB 中 Mn 未与 HEXXH 配位（常见：Mn 坐标被移到 DNA 磷酸附近、或清理 PDB 时丢失了 Mn → 从 holo 参考结构叠加重放 Mn）。
- **铁律**：Mn 三齿配位（PprI: H92 NE2 / H96 NE2 / E123 OE1+OE2）是催化态起点；几何不对后面 DISANG 会把体系拉变形，prod 首帧 RESTRAINT 能量会 >100（§5）。

## 3. frcmod（12-6 非键参数）

- 仓库自带 `ff/mn_cm12-6.frcmod`（CM 12-6 Mn²⁺，Li et al. *JCTC* 2013, 9, 2733; DOI 10.1021/ct400146w）。来源与备选参数见 `ff/README_FF.md`。
- 定位逻辑：`METAL=yes` 时 `FRCMOD` 留空 → 自动用仓库 `ff/mn_cm12-6.frcmod`；文件缺失即 `die`。
- stage1 tleap 中 `loadamberparams` 加载；拓扑校验（parmed）会确认 **Mn 电荷 = +2**，否则退出。

## 4. masks.env / DISANG（stage2 产物，自动生成）

- `make_restraints.py` 读 parm7/rst7 自动识别蛋白/DNA/Mn 区间与配位原子，输出到 `topol/<SYSTEM>/masks.env`：
  - `PROT_RANGE` / `DNA_RANGE`（残基区间，如 `1-254` / `255-278`）
  - `MN_RES`（如 `279`；`-` = 未识别到 Mn）
  - `COORD1..4`（如 `95:NE2` = 残基95 NE2 原子，即 PprI 的 H92——**注意 parm7 残基号 = 蛋白残基 + 1 的位移关系**，DNA 从 255 起）
- 检查项：MN_RES 非 `-`；COORD1..4 非空；若配位原子数与预期不符（PprI 应 4 个：H92/H96/E123×2 羧基氧），看 make_restraints 日志确认探测半径（默认 3.5 Å）。
- **无金属体系**：`METAL=no` → stage2 写空 disang.txt 占位，阶段 3 模板不启用 nmropt。

## 5. prod 首帧 RESTRAINT 能量闸门（stage3）

- prod 模板只保留 Mn 配位 flat-bottom DISANG（**2.1–2.5 Å，rk=10 kcal/mol·Å²**），骨架 **ntr=0 自由**（prod 铁律，见 README）。
- stage3 读 prod.out 首帧 RESTRAINT 能量：**健康 <6 kcal/mol**；>100 → 报错退出并提示本清单。
- >100 常见原因：起始 DNA（或水）贴近 Mn 造成强 repulsion、配位原子选错（DISANG 约束到非配位原子）、起始几何已崩（§2）。

## 6. 运行期健康信号（模拟结束后核对）

| 量 | 预期 | 来源 |
|---|---|---|
| Mn²⁺ 配位数 | ~6（三齿 + 水/磷酸），**>7 说明 12-6 过结合** | cpptraj `an5_rdf.in` + 配位距离时序 |
| 配位距离时序 (fig08) | H92/H96 严格保持 ~2.0–2.4 Å；E123 双齿可微波动 | `results/<SYSTEM>/mn_H71/mn_H75/mn_E102a/mn_E102b.dat` |
| Mn RDF 第一壳层 (fig30) | ~1.95–2.1 Å 主峰（12-6 模型） | `an5_rdf` |
| 蛋白 RMSD (fig01) | 1–3 Å 平台（骨架自由后正常涨落） | `rmsd_prot.dat` |

## 7. 诚实边界（MM-PBSA / 绝对能量）

- 12-6 非键模型**低估**二价阳离子-配体相互作用（HFE 偏差 ~ −24 kcal/mol）；本流水线把 Mn 用在「口袋内 + DISANG 锚定」的相对比较场景，误差可接受。
- **MM-PBSA 剥 Mn**（stage5 已自动处理）：sander GB 不识别 Mn 原子类型 → rec/lig/cp 拓扑与轨迹都剥 Mn，Mn 不参与 ΔG。
- 判断设计好坏看**相对差**（同一蛋白不同 DNA 的 ΔΔG / 结合模式），不看 Mn 绝对能量。
