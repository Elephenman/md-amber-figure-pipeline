# 论文 Methods 模板（英文，可直接改写投稿）

> 对应 `md-amber-figure-pipeline` 默认参数与协议（PprI·ssDNA·Mn²⁺ 案例）。
> 含三段：① AMBER 选型与体系构建辩护 ② 金属 12-6 参数化与配位约束 ③ 生产 MD 与后处理。
> 替换 `[...]` 内占位符（蛋白名 / PDB / 时长 / 重复数 / 帧数）后即可用。
> 中文批注用 `【】` 标出，投稿前删除。

---

## 1. System building and force-field choice

All molecular dynamics (MD) simulations were performed with the **AMBER 22** package [Case et al.]. We deliberately chose AMBER over GROMACS for three reasons: (i) the AMBER `pmemd` engine natively supports **12-6 Lennard-Jones non-bonded parameters** for divalent cations, which is not the case for the 12-6-4 model without patched support; (ii) the built-in `cpptraj` and `MMPBSA.py` tools provide an integrated analysis path for per-residue decomposition of binding free energies; and (iii) the ff19SB/OL15 parameter combination is the current recommended standard for protein–nucleic-acid complexes in AMBER.

【选型辩护三理由：12-6 原生支持 / cpptraj+MMPBSA 一体化 / ff19SB+OL15 组合标准】

The initial structure of the [protein–ssDNA–Mn²⁺] complex was taken from [PDB: 8SLN] ... Missing loops were modeled [with ...]; hydrogen atoms were added with `tleap`. The protein and the single-stranded DNA were described by the **ff19SB** [Tian et al.] and **OL15** [Zgarbová et al.] force fields, respectively, and the system was solvated in a truncated octahedron TIP3P water box [Jorgensen et al.] with a 12 Å buffer. Sodium and chloride counterions (12-6-4 LJ parameters, `ions1lm_1264_tip3p`) were added to neutralize the system to a target ionic strength of [0.15] M.

【力场/溶剂/离子：ff19SB+OL15+TIP3P+12-6-4 Na+/Cl-，12 Å 截断】

## 2. Mn²⁺ parameterization and coordination restraints

Mn²⁺ was modeled with the **CM 12-6 parameters** of Li et al. [*J. Chem. Theory Comput.* 2013, 9, 2733–2748], with R*/2 = 1.407 Å, ε = 0.016867 kcal/mol and charge +2, loaded through a `frcmod` file. A 12-6 (rather than 12-6-4) model was adopted because production runs use `pmemd.cuda`, which has no native 1/r⁴ term; within the metal-binding pocket, where the ion is additionally restrained, the known underestimation of 12-6 models toward divalent-cation–ligand interactions (hydration free-energy bias ~ −24 kcal/mol) is tolerable because our readouts are *relative* (same protein, different DNA).

【12-6 参数来源与诚实边界：口袋内+约束场景、只读相对差】

Histidine residues coordinating the metal were assigned the correct protonation states by an automated protocol: HIS residues within 3.0 Å of Mn²⁺ were rewritten as HID (ε-protonated/δ-deprotonated) or HIE depending on the coordinating atom. To preserve the catalytic geometry [His92/His96/Glu123 tridentate coordination; GAF/peptidase dyad in the PprI case], the Mn²⁺–ligand distances were restrained with **flat-bottom harmonic restraints of 2.1–2.5 Å (force constant 10 kcal·mol⁻¹·Å⁻²)** applied during equilibration and production. During production, no other positional restraints were applied (ntr=0), because even weak backbone restraints (wt=1.0) were found to pin the protein RMSD at ~0.39 Å and produce pathologically limited sampling.

【配位 His→HID/HIE 自动改写；flat-bottom DISANG 2.1–2.5 Å rk=10；prod ntr=0 自由 —— prod 铁律】

## 3. Production MD protocol

Each system was minimized in two steps — (i) 3,000 cycles with backbone restraints of 25 kcal·mol⁻¹·Å⁻², and (ii) 5,000 unrestrained cycles — followed by heating from 0 to 300 K over 100 ps in the NVT ensemble with backbone restraints of 10 kcal·mol⁻¹·Å⁻², and a 500 ps NPT equilibration at 300 K with backbone restraints of 2 kcal·mol⁻¹·Å⁻². Production runs of **[150] ns** per replicate (frames saved every 100 ps → 1,500 frames) were carried out at 300 K in the NPT ensemble with the Langevin thermostat and Monte Carlo barostat, a 8 Å non-bonded cutoff, and PME electrostatics. **[N] independent replicates** were performed starting from ...; results are reported as mean ± SD over replicates where appropriate.

【min1(wt25,3000cyc)→min2(5000cyc 无约束)→heat(NVT 100ps wt10)→equil(NPT 500ps wt2)→prod(NPT ntr=0, 100ps/帧)。写实际值替换 150/1/N】

Trajectory analyses (RMSD, RMSF, radius of gyration, solvent-accessible surface area, PCA, DSSP, inter-domain dynamic cross-correlation matrices, and Mn²⁺ coordination distances / radial distribution functions) were performed with `cpptraj`. Binding free energies were estimated by the **MM-PBSA/GBSA** approach (`MMPBSA.py`, GB model igb=5 with salt concentration 0.15 M for per-residue decomposition; PB with ionic strength 0.15 M for total energies), sampling every 10th frame of the production trajectory. Because the sander GB engine does not recognize the Mn atom type, Mn²⁺ was stripped from receptor, ligand and complex topologies and trajectories before all MM-PBSA calculations. Figures and galleries were rendered with in-house Python scripts (numpy/scipy/matplotlib) [repository: github.com/Elephenman/md-amber-figure-pipeline].

【后处理：cpptraj 全量量 + MM-PBSA(GB igb=5 per-residue + PB 总量, interval=10) + 剥 Mn 铁律】

---

## 引用建议（投稿时按期刊规范补全）

- AMBER 22: Case, D.A. et al. *J. Chem. Inf. Model.* 2023 (AmberTools 23 / AMBER 22 cite)
- ff19SB: Tian, C. et al. *J. Chem. Theory Comput.* 2020, 16, 528–552
- OL15: Zgarbová, M. et al. *J. Chem. Theory Comput.* 2015, 11, 5723–5736
- TIP3P: Jorgensen, W.L. et al. *J. Chem. Phys.* 1983, 79, 926–935
- CM 12-6 Mn²⁺: Li, P.; Roberts, B.P.; Chakravorty, D.K.; Merz, K.M. Jr. *J. Chem. Theory Comput.* 2013, 9, 2733–2748
- MMPBSA.py: Miller, B.R. III et al. *J. Chem. Theory Comput.* 2012, 8, 3314–3321

> 数值摘要（PprI WT·S1 案例，供比对量级）：蛋白 RMSD 1.93±0.36 Å；Mn 三齿配位 H92/H96 全程保持；界面 H 键 R253–T24 occupancy 76.7%；MM-PBSA ΣΔG −126.2 kcal/mol，R207 单残基贡献 −13.0 kcal/mol；MD RMSF 与晶体 B-factor Pearson r = 0.62。
