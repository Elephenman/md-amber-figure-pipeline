#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fix_his.py — 通用金属配位 His 质子化态修正 (HIS → HID/HIE)

原理: Mn²⁺ 等金属常配位 His 咪唑 N（NE2=εN 或 ND1=δN）。
  配位 N 必须去质子化(不带 H)才能结合金属:
    - Mn 配位 NE2 → εN 空出、δN 带 H  → 残名改 HID
    - Mn 配位 ND1 → δN 空出、εN 带 H  → 残名改 HIE
  （ff19SB leaprc 只认 HID/HIE/HIP，旧式 HSD/HSE/HSP 一并转写）
  不改则 tleap 默认 HIS→HIE (εN 带 H)，与金属 ~2A 处 clash。

用法: python3 fix_his.py <input.pdb> [output.pdb] [--thresh 3.0] [--chain X]
  - 未给 output.pdb → 原地覆写（先备份 <input>.bak）
  - 找不到 Mn 或无配位 His → 不改动，退出码 0
输出: 打印每个改写残基行号; 总改数
"""
import os
import sys

THRESH = 3.0  # Mn–N 配位判距 (Å)


def load_coords(path):
    """解析 ATOM/HETATM -> [(line_idx, name, resn, chain, resi, elem, x,y,z)]"""
    rows = []
    with open(path, encoding="utf-8", errors="ignore") as f:
        for i, line in enumerate(f):
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
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            elem = line[76:78].strip()
            rows.append((i, name, resn, chain, resi, elem, x, y, z))
    return rows


def main():
    if len(sys.argv) < 2:
        sys.stderr.write("用法: python3 fix_his.py <input.pdb> [output.pdb] [--thresh 3.0] [--chain X]\n")
        return 2
    inp = sys.argv[1]
    out = None
    chain_only = None
    for i, a in enumerate(sys.argv[2:], start=2):
        if a == "--thresh":
            global THRESH
            THRESH = float(sys.argv[i + 1])
        elif a == "--chain":
            chain_only = sys.argv[i + 1]
        elif not a.startswith("--") and out is None:
            out = a
    if not os.path.isfile(inp):
        sys.stderr.write(f"输入不存在: {inp}\n")
        return 1
    rows = load_coords(inp)
    mns = [r for r in rows if r[5] == "MN"]
    if not mns:
        print(f"[info] {inp}: 无 Mn 原子 — 跳过 HID/HIE 改写")
        return 0

    # 组残基 → 取 NE2/ND1 坐标
    his = {}
    for r in rows:
        if r[2] in ("HIS", "HSD", "HSE", "HSP") and r[0] >= 0 and r[4] > 0:
            his.setdefault((r[3], r[4]), []).append(r)
    if not his:
        print(f"[info] {inp}: 有 Mn 但无 His — 跳过")
        return 0

    # 决定每个 His 改什么: 统计 Mn 距 NE2/ND1 最近者
    to_fix = {}  # (chain, resi) -> "HID"/"HIE"
    for (chain, resi), rs in his.items():
        if chain_only and chain != chain_only:
            continue
        nd1 = next((r for r in rs if r[1] == "ND1"), None)
        ne2 = next((r for r in rs if r[1] == "NE2"), None)
        for mn in mns:
            best = None  # (dist, target_n)
            for cand, tag in ((nd1, "HIE"), (ne2, "HID")):
                if cand is None:
                    continue
                d = ((cand[6] - mn[6]) ** 2 + (cand[7] - mn[7]) ** 2 + (cand[8] - mn[8]) ** 2) ** 0.5
                if best is None or d < best[0]:
                    best = (d, tag, cand[1])
            if best and best[0] < THRESH:
                # 配位 NE2 -> εN 去质子 = HID; 配位 ND1 -> δN 去质子 = HIE
                to_fix[(chain, resi)] = best[1]
                break

    if not to_fix:
        print(f"[info] {inp}: Mn 与 His 距离均 >{THRESH:.1f} Å — 不改写")
        return 0

    # 写回（原地先备份）
    with open(inp, encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()
    n_changed = 0
    for i, line in enumerate(lines):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resn = line[17:20]
        if resn not in ("HIS", "HSD", "HSE", "HSP"):
            continue
        chain = line[21].strip() or "A"
        try:
            resi = int(line[22:26])
        except ValueError:
            continue
        newn = to_fix.get((chain, resi))
        if newn and line[17:20] != newn:
            lines[i] = line[:17] + newn + line[20:]
            n_changed += 1

    dst = out or inp
    if out is None and os.path.exists(inp):
        try:
            os.replace(inp, inp + ".bak")
        except OSError:
            pass
    with open(dst, "w", encoding="utf-8", newline="\n") as f:
        f.writelines(lines)
    for (chain, resi), tag in sorted(to_fix.items()):
        print(f"  {chain}{resi}: HIS → {tag}  (Mn 配位)")
    print(f"[OK] {inp} → {dst}: 改写 {n_changed} 条 ATOM 记录")
    return 0


if __name__ == "__main__":
    sys.exit(main())
