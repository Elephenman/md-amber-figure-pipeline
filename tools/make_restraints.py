#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_restraints.py — 由 parm7/rst7 生成：
  1) AMBER 约束 mask（蛋白/DNA 骨架 + Mn）
  2) DISANG 距离约束（Mn²⁺ → HEXXH 配位原子，flat-bottom）
用法: python3 make_restraints.py <parm7> <rst7> <outdir>
输出: outdir/masks.env  (HEAVY_BB_MASK / BB_MASK)
      outdir/disang.txt (DISANG 距离约束)
"""
import os
import sys
import parmed as pmd

# 2026-09-01: 输入 PDB 的 H71/H75 已 HID 化(HSD) -> prmtop 残名 HID，须纳入蛋白集合，
# 否则蛋白残基区间漏匹配、配位原子选择失败（原仅 "HIS" 匹配不到 HID）。
PROT_AA = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","HID","HIE","HIP","ILE",
           "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
DNA_RES = {"DA","DC","DG","DT","DA5","DT5","DA3","DT3","DC5","DG5","DC3","DG3",
            "DU","DU5","DU3"}  # 含 5'/3' 末端变体（parmed mask :DA 不匹配 DA5/DT3）
PRIMES = {"P","OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'"}  # DNA 骨架原子
MN_THRES = 3.5  # 配位原子判距（Å）


def main():
    parm7, rst7, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(outdir, exist_ok=True)
    parm = pmd.load_file(parm7, rst7)

    prot_ids, dna_ids, mn_ids = [], [], []
    mn_atoms = []
    for res in parm.residues:
        name = res.name
        if name == "MN":
            mn_ids.append(res.idx)
            mn_atoms.extend(res.atoms)
        elif name in PROT_AA:
            prot_ids.append(res.idx + 1)   # AMBER mask 残基编号 1-based
        elif name in DNA_RES:
            dna_ids.append(res.idx + 1)   # AMBER mask 残基编号 1-based
        else:  # 溶剂/离子
            pass

    # 2026-09-01 REFINE_V2: 支持三种体系类型（complex / apo / freedna）
    #   freedna 无 Mn -> 跳过 DISANG；apo 无 DNA、freedna 无蛋白 -> 区间按实际。
    if not mn_atoms:
        print("[WARN] 无 Mn 残基（freedna/apo 金属-free）— 跳过 DISANG 约束", file=sys.stderr)
        mn_idx, mn_res = None, ""
    else:
        mn = mn_atoms[0]
        mn_idx = mn.idx + 1          # 1-based 原子序号（用于 DISANG iat）
        mn_res = mn.residue.idx + 1  # 1-based 残基序号（用于 MMPBSA mask :MN_RES）

    def mask_for(ids, atom_sel):
        return "(" + " | ".join(f"({r}:{atom_sel})" for r in ids) + ")"
        # 实际用区间更紧凑:
    # 残基区间生成：直接使用连续区间（parm7 中蛋白 1..254、DNA 255..278 连续分块）。
    # 不要用碎片化区间（如 :1-24,:26-50）：AMBER 旧版 atommask 解析器对重复冒号
    # 的逗号列表支持差（"unknown symbol ::"）；连续单区间最稳。
    # 若 PROT_AA 漏匹配（HID/HIE 等质子化变体）导致区间碎片化，这里强制用连续区间。
    prot_rng, dna_rng = "", ""
    if prot_ids:
        prot_lo, prot_hi = min(prot_ids), max(prot_ids)
        prot_rng = f":{prot_lo}-{prot_hi}"
    if dna_ids:
        dna_lo, dna_hi = min(dna_ids), max(dna_ids)
        dna_rng = f":{dna_lo}-{dna_hi}"

    # 蛋白骨架 = CA,C,N,O；DNA 骨架 = 磷酸/糖主链；Mn 单独（按存在性拼接）
    bb_parts, bb2_parts = [], []
    if prot_rng:
        bb_parts.append(f"(({prot_rng}) & @CA,C,N,O)")
        bb2_parts.append(f"(({prot_rng}) & @CA,C,N,O)")
    if dna_rng:
        bb_parts.append(f"(({dna_rng}) & @P,OP1,OP2,O5',C5',C4',O4',C3',O3')")
        bb2_parts.append(f"(({dna_rng}) & @P,OP1,OP2,O5',C5',C4',O4',C3',O3')")
    if mn_res:
        bb_parts.append(":MN")
    HEAVY_BB = " | ".join(bb_parts) if bb_parts else ":WAT"
    BB = " | ".join(bb2_parts) if bb2_parts else ":WAT"

    # 配位原子：Mn 3.5A 内、蛋白侧链（HIS/HID/HIE 的 NE2/ND1、GLU 的 OE1/OE2）
    coords = []
    if mn_atoms:
        mn = mn_atoms[0]
        for a in parm.atoms:
            if a.idx == mn.idx or a.residue.name in ("WAT", "Na+", "Cl-"):
                continue
            if a.residue.name in ("HIS", "HID", "HIE", "HIP") and a.name in ("NE2", "ND1"):
                coords.append((a, "HIS"))
            elif a.residue.name == "GLU" and a.name in ("OE1", "OE2"):
                coords.append((a, "GLU"))
        # 用最小化后坐标算距离
        mn_xyz = mn.xx, mn.xy, mn.xz
        chosen = []
        for a, kind in coords:
            d = ((a.xx - mn_xyz[0]) ** 2 + (a.xy - mn_xyz[1]) ** 2 + (a.xz - mn_xyz[2]) ** 2) ** 0.5
            if d < MN_THRES:
                chosen.append((a, d, kind))
        chosen.sort(key=lambda t: t[1])
        if not chosen:
            print("[WARN] 未找到 3.5A 内 Mn 配位蛋白原子，将只保留 Mn 位置约束", file=sys.stderr)
    else:
        chosen = []

    # DISANG：Mn–配位原子 flat-bottom（平区 2.1–2.5A）
    with open(os.path.join(outdir, "disang.txt"), "w") as f:
        f.write("Flat-bottom distance restraints: Mn2+ to HEXXH coordination atoms\n")
        for a, d, kind in chosen[:4]:
            idx2 = a.idx + 1
            f.write("&rst\n")
            f.write(f"  iat={mn_idx},{idx2},\n")
            f.write("  r1=1.9, r2=2.1, r3=2.5, r4=2.7, rk2=10.0, rk3=10.0,\n")
            f.write("/\n")
        print(f"[disang] 写入 {len(chosen[:4])} 条 Mn 配位约束"
              + (f"（{mn.residue.name}{mn.residue.idx+1}）" if mn_atoms else "（无 Mn）"))

    # 输出 masks（PROT_RANGE/DNA_RANGE/MN_RES 按实际；空值用 "-" 占位）
    # COORD1..4: Mn 3.5 A 内配位原子 "AMBER残基号:原子名"（按距离升序，最多 4）
    with open(os.path.join(outdir, "masks.env"), "w") as f:
        f.write(f"export HEAVY_BB_MASK=\"{HEAVY_BB}\"\n")
        f.write(f"export BB_MASK=\"{BB}\"\n")
        f.write(f"export PROT_RANGE=\"{prot_rng[1:] if prot_rng else '-'}\"\n")
        f.write(f"export DNA_RANGE=\"{dna_rng[1:] if dna_rng else '-'}\"\n")
        f.write(f"export MN_RES=\"{mn_res if mn_res else '-'}\"\n")
        for i in range(4):
            if i < len(chosen):
                a, _, _ = chosen[i]
                f.write(f"export COORD{i+1}=\"{a.residue.idx + 1}:{a.name}\"\n")
            else:
                f.write(f"export COORD{i+1}=\"-\"\n")
    print("HEAVY_BB_MASK:", HEAVY_BB)
    print("BB_MASK:", BB)
    print(f"PROT_RANGE={prot_rng[1:] if prot_rng else '-'} DNA_RANGE={dna_rng[1:] if dna_rng else '-'} MN_RES={mn_res or '-'}")
    print("COORD1..4:", [(chosen[i][0].residue.idx + 1, chosen[i][0].name) if i < len(chosen) else "-" for i in range(4)])


if __name__ == "__main__":
    main()
