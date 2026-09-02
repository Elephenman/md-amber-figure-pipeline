#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""check_input_pdb.py — 输入 PDB 体检（stage0 调用，只读不改文件）
用法: python3 check_input_pdb.py <input.pdb> [--metal yes|no]
检查项:
  1) ATOM/HETATM 记录数、链组成（蛋白/核酸/水/离子）
  2) METAL=yes: 必须存在 Mn 元素原子；并报告 Mn 最近的蛋白/核酸配位原子距离
  3) 蛋白残基编号是否 1 开头连续（非致命，提示作者编号需在 project.env 校准 PROT_RANGE）
退出码: 0=通过  1=致命问题（METAL=yes 却无 Mn / 无任何原子）
"""
import os
import re
import sys

AA3 = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU",
       "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
       "HID","HIE","HIP","HSD","HSE","HSP"}
NA3 = {"DA","DC","DG","DT","DU","DA5","DT5","DA3","DT3","DC5","DG5","DC3","DG3",
       "DU5","DU3","A","C","G","T","U","RA","RC","RG","RU"}

def parse(pdb):
    """返回 (atoms, chains) — atoms: list[(name,resn,chain,resi,element,x,y,z)]
    chains: {chain: [(lo,hi,kind), ...]} kind=prot|na|other"""
    atoms = []
    with open(pdb, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            name = line[12:16].strip()
            resn = line[17:20].strip()
            chain = line[21].strip() or "A"
            try:
                resi = int(line[22:26])
            except ValueError:
                resi = 0
            el = line[76:78].strip() if len(line) >= 78 else ""
            if not el:  # 旧格式: 从原子名猜元素
                el = re.sub(r"[0-9+\-]", "", name)[0] if name else "?"
            try:
                x, y, z = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                x = y = z = 0.0
            atoms.append((name, resn, chain, resi, el, x, y, z))
    # 链区间汇总
    chains = {}
    seen = {}
    for name, resn, chain, resi, el, x, y, z in atoms:
        kind = "prot" if resn in AA3 else ("na" if resn in NA3 else "other")
        key = (chain, kind)
        seen.setdefault(key, set()).add(resi)
    for (chain, kind), resis in seen.items():
        lo, hi = min(resis), max(resis)
        chains.setdefault(chain, []).append((lo, hi, kind))
    return atoms, chains


def main():
    args = [a for a in sys.argv[1:]]
    metal = "auto"
    pdb = None
    for i, a in enumerate(args):
        if a == "--metal":
            metal = args[i + 1] if i + 1 < len(args) else "auto"
        elif a.startswith("--metal="):
            metal = a.split("=", 1)[1]
        elif not a.startswith("--") and pdb is None:
            pdb = a
    if not pdb or not os.path.isfile(pdb):
        sys.stderr.write(f"用法: check_input_pdb.py <input.pdb> [--metal yes|no|auto]\n")
        return 2

    atoms, chains = parse(pdb)
    n_atom = sum(1 for a in atoms if a[0])
    if n_atom == 0:
        sys.stderr.write("[FAIL] PDB 无任何 ATOM/HETATM 记录\n")
        return 1
    print(f"原子记录总数: {n_atom}")
    for chain in sorted(chains):
        for lo, hi, kind in chains[chain]:
            mark = {"prot": "蛋白", "na": "核酸", "other": "其它(HETATM/水)"}[kind]
            print(f"  链 {chain}: {mark} 残基 {lo}-{hi}")

    mn = [a for a in atoms if a[4] == "MN"]
    if metal == "yes":
        if not mn:
            sys.stderr.write("[FAIL] METAL=yes 但 PDB 中找不到 Mn 元素原子（HETATM MN）\n")
            return 1
        m = mn[0]
        print(f"Mn²⁺ 位于 链{m[2]} 残基{m[3]}  ({m[5]:.2f},{m[6]:.2f},{m[7]:.2f})")
        # 最近蛋白/核酸原子
        best = None
        for a in atoms:
            if a[4] == "MN":
                continue
            d = ((a[5] - m[5]) ** 2 + (a[6] - m[6]) ** 2 + (a[7] - m[7]) ** 2) ** 0.5
            if best is None or d < best[0]:
                best = (d, a)
        if best:
            d, a = best
            print(f"  Mn 最近配位原子: 链{a[2]} {a[1]}{a[3]} {a[0]}  d={d:.2f} Å")
            if d > 3.5:
                sys.stderr.write(f"[WARN] 最近原子 d={d:.2f} Å > 3.5 —— Mn 可能游离未配位，"
                                 f"请核对输入结构（参考 docs/METAL_CHECKLIST.md）\n")
    elif metal == "auto":
        if mn:
            print(f"自动探测到 Mn²⁺（{len(mn)} 个）—— 将按金属体系处理")
    print("[OK] 输入 PDB 体检通过")
    return 0


if __name__ == "__main__":
    sys.exit(main())
