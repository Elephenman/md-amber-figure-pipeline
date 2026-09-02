#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
figures_gallery.py — 通用 HTML 画廊生成器（md-amber-figure-pipeline）
扫描 <MDEASY_RES>/figures_pretty/fig*.png 输出 index.html（深色科研风单文件）。

用法:
    python figures_gallery.py                 # 读 env MDEASY_RES（stage6 已导出）
    MDEASY_RES=.../results/<SYSTEM> python figures_gallery.py

分组按 fig 编号段: 01-08 基础 / 09-22 深入机理 / 23-30 论文级(图集 III)
无需手工维护条目: 标题/分组从文件名自动派生。
"""
import os, re, html, glob

def _norm_msys(p):
    """Git Bash(MSYS) 路径 /a/... → A:/...（Windows 原生 python 可读）"""
    m = re.match(r"^/([a-zA-Z])/(.*)$", p)
    return f"{m.group(1).upper()}:/{m.group(2)}" if m else p

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = _norm_msys(os.environ.get("MDEASY_RES") or
                    os.path.join(HERE, "results", "WT__S1"))
FIGDIR = os.path.join(FIGDIR, "figures_pretty")

GROUPS = [("基础", "01"), ("深入", "09"), ("论文级", "23")]


def group_of(fig):
    n = int(fig[:2])
    for gname, start in reversed(GROUPS):
        if n >= int(start):
            return gname
    return "其他"


def title_of(fname):
    t = re.sub(r"^fig\d+_|\.png$", "", fname).replace("_", " ").strip()
    return t or fname


def esc(s):
    return html.escape(str(s))


def main():
    files = sorted(glob.glob(os.path.join(FIGDIR, "fig*.png")))
    if not files:
        print(f"[skip] {FIGDIR} 无 fig*.png")
        return 0
    cards = []
    for p in files:
        fn = os.path.basename(p)
        m = re.match(r"fig(\d+)_", fn)
        fig = m.group(1) if m else "?"
        rel = os.path.relpath(p, FIGDIR).replace("\\", "/")
        cards.append((group_of(fig), fig, rel, title_of(fn)))
    cards.sort(key=lambda c: (c[0] != "基础", int(c[2][3:5]) if c[2][3:5].isdigit() else 99))

    grouped = {}
    for g, fig, rel, title in cards:
        grouped.setdefault(g, []).append((fig, rel, title))

    card_html = []
    for gname in [g[0] for g in GROUPS] + [k for k in grouped if k not in [g[0] for g in GROUPS]]:
        if gname not in grouped:
            continue
        items = "".join(
            f'<div class="card"><img loading="lazy" src="{esc(rel)}" alt="fig{fig}">'
            f'<div class="cap"><span class="tag">fig{fig}</span>{esc(title)}</div></div>'
            for fig, rel, title in grouped[gname])
        card_html.append(f'<section><h2>{esc(gname)}</h2><div class="grid">{items}</div></section>')

    doc = f"""<!DOCTYPE html>
<html lang="zh"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>MD 图集 — {esc(os.path.basename(os.path.dirname(FIGDIR)))}</title>
<style>
:root {{ --bg:#0E1117; --panel:#161B22; --border:#2A3140; --fg:#E6E8EC;
        --dim:#8B949E; --acc:#58A6FF; }}
* {{ box-sizing:border-box; margin:0; padding:0; }}
body {{ background:var(--bg); color:var(--fg);
       font-family:"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif; }}
header {{ padding:36px 48px 22px; border-bottom:1px solid var(--border);
         background:linear-gradient(180deg,#10151D,#0E1117); }}
header h1 {{ font-size:22px; font-weight:700; }}
header small {{ color:var(--dim); font-weight:400; margin-left:10px; }}
.meta {{ color:var(--dim); font-size:12.5px; margin-top:8px; line-height:1.8; }}
main {{ padding:26px 48px 60px; }}
section {{ margin-top:28px; }}
h2 {{ font-size:15px; color:var(--acc); letter-spacing:1px; margin-bottom:12px;
      border-left:3px solid var(--acc); padding-left:10px; }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fill,minmax(300px,1fr));
        gap:16px; }}
.card {{ background:var(--panel); border:1px solid var(--border); border-radius:10px;
        overflow:hidden; transition:transform .15s, border-color .15s; }}
.card:hover {{ transform:translateY(-2px); border-color:var(--acc); }}
.card img {{ width:100%; display:block; }}
.cap {{ padding:9px 12px; font-size:12.5px; color:var(--dim); }}
.tag {{ color:var(--acc); font-weight:700; margin-right:8px; }}
footer {{ padding:20px 48px 40px; color:var(--dim); font-size:12px;
         border-top:1px solid var(--border); line-height:1.9; }}
</style></head><body>
<header>
  <h1>Amber MD 图集<small> · {len(cards)} 张</small></h1>
  <div class="meta">体系: {esc(os.path.basename(os.path.dirname(FIGDIR)))}<br>
  流水线: <code>md-amber-figure-pipeline</code> — stage4 cpptraj 采集 → stage5 MM-PBSA(可选) → stage6 一键出图<br>
  工具链: AmberTools cpptraj · MMPBSA.py (GB) · numpy/scipy/matplotlib</div>
</header>
<main>{''.join(card_html)}</main>
<footer>由 scripts/figures_gallery.py 自动生成（扫描 figures_pretty/fig*.png）。</footer>
</body></html>"""
    out = os.path.join(FIGDIR, "index.html")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(doc)
    print(f"[written] {out}  ({len(doc):,} bytes, {len(cards)} figs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
