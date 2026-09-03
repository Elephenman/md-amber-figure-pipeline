# -*- coding: utf-8 -*-
"""
make_standalone_gallery.py — 把 figures_pretty/ 的 index.html + N 张 PNG 打包成
**单个自包含 HTML**（零外部依赖、零外链、可离线双击打开、可直接分享）。

设计：
  1) 网格缩略图  = 轻量 JPEG (默认宽 800px, q82) 内嵌 base64 → 页面秒开
  2) 点击卡片     = 灯箱加载 **原图无损 PNG** 内嵌 base64 → 出版级放大
  3) 支持 键盘 ←/→ 翻页、Esc 关闭、点击背景关闭、单图下载、文字即时筛选

用法：
  python make_standalone_gallery.py                       # 默认读 ./results/WT__S1/figures_pretty
  python make_standalone_gallery.py --dir <figures_dir> --out <file.html>
  python make_standalone_gallery.py --thumb-width 1000 --thumb-quality 88
  python make_standalone_gallery.py --full jpeg --full-quality 92   # 瘦身版(有损)

输出默认: <figures_dir>/../<CASE>_figures_standalone.html
"""
import argparse, base64, io, os, re, sys, datetime

try:
    from PIL import Image
except ImportError:
    sys.exit("[err] 需要 Pillow: pip install pillow")


# ---------------------------------------------------------------- utils
def b64_png(path):
    with open(path, "rb") as f:
        return "data:image/png;base64," + base64.b64encode(f.read()).decode("ascii")


def b64_thumb(path, width, quality):
    """RGBA -> 白底 RGB -> 等比缩放到 width -> JPEG base64."""
    im = Image.open(path)
    if im.mode in ("RGBA", "LA", "P"):
        bg = Image.new("RGB", im.size, (255, 255, 255))
        im = im.convert("RGBA")
        bg.paste(im, mask=im.split()[-1])
        im = bg
    else:
        im = im.convert("RGB")
    if im.width > width:
        h = round(im.height * width / im.width)
        im = im.resize((width, h), Image.LANCZOS)
    buf = io.BytesIO()
    im.save(buf, "JPEG", quality=quality, optimize=True, progressive=True)
    return "data:image/jpeg;base64," + base64.b64encode(buf.getvalue()).decode("ascii")


def b64_full(path, kind, quality):
    if kind == "png":
        return b64_png(path)
    im = Image.open(path)
    if im.mode in ("RGBA", "LA", "P"):
        bg = Image.new("RGB", im.size, (255, 255, 255))
        im = im.convert("RGBA")
        bg.paste(im, mask=im.split()[-1])
        im = bg
    else:
        im = im.convert("RGB")
    buf = io.BytesIO()
    im.save(buf, "JPEG", quality=quality, optimize=True, progressive=True, subsampling=0)
    return "data:image/jpeg;base64," + base64.b64encode(buf.getvalue()).decode("ascii")


EXTRA_CSS = """
  /* ---- standalone gallery additions ---- */
  .toolbar { padding: 12px 48px; background: #101821; border-bottom: 1px solid var(--border);
             display: flex; flex-wrap: wrap; gap: 10px; align-items: center; }
  .toolbar input { flex: 1 1 260px; min-width: 200px; padding: 8px 12px; font-size: 13px;
                   border-radius: 4px; border: 1px solid #2A3648; background: #0E141C;
                   color: var(--ink); outline: none; }
  .toolbar input:focus { border-color: var(--accent); }
  .toolbar .hint { font-size: 12px; color: var(--dim); }
  .toolbar .count { font-size: 12.5px; color: #7AB7E8; font-family: ui-monospace,Consolas,monospace; }
  .card { cursor: zoom-in; position: relative; }
  .card.hidden { display: none; }
  .card .zoom { position: absolute; top: 12px; right: 14px; font-size: 11px; color: var(--dim);
                background: rgba(14,20,28,.85); padding: 2px 7px; border-radius: 3px;
                opacity: 0; transition: opacity .15s; }
  .card:hover .zoom { opacity: 1; }
  .grp-title.hidden { display: none; }

  #lb { position: fixed; inset: 0; background: rgba(6,9,13,.94); z-index: 999;
        display: none; flex-direction: column; }
  #lb.on { display: flex; }
  #lb .bar { display: flex; align-items: center; gap: 14px; padding: 12px 22px;
             background: #131A24; border-bottom: 1px solid var(--border); flex: 0 0 auto; }
  #lb .bar .no { font-weight: 700; font-size: 15px; color: #FFF;
                 font-family: ui-monospace,Consolas,monospace; }
  #lb .bar .ttl { font-size: 14px; color: var(--ink); }
  #lb .bar .sp { flex: 1; }
  #lb .bar button { background: #1E2837; color: var(--ink); border: 1px solid #2A3648;
                    border-radius: 4px; padding: 6px 13px; font-size: 13px; cursor: pointer; }
  #lb .bar button:hover { background: #2A3648; }
  #lb .stage { flex: 1 1 auto; overflow: auto; display: flex; align-items: center;
               justify-content: center; padding: 18px; }
  #lb .stage img { max-width: 100%; max-height: 100%; background: #FFF;
                   box-shadow: 0 8px 40px rgba(0,0,0,.6); }
  #lb .foot { flex: 0 0 auto; padding: 10px 22px 16px; text-align: center;
              font-size: 12.2px; color: var(--dim); background: #131A24;
              border-top: 1px solid var(--border); line-height: 1.6; }
  #lb .foot b { color: #E8623C; font-weight: 600; }
  @media (max-width: 720px) {
    header, .toc, main, .toolbar { padding-left: 18px; padding-right: 18px; }
    .grid { grid-template-columns: 1fr; }
  }
"""

EXTRA_JS = r"""
(function(){
  var FULL = window.__FULL__ || [];
  var cards = [].slice.call(document.querySelectorAll('.card'));
  var meta  = cards.map(function(c){
    var t = c.querySelector('.title'), k = c.querySelector('.kv');
    return { idx: +c.dataset.idx,
             no:  (c.querySelector('.fig-no')||{}).textContent || '',
             ttl: t ? t.textContent : '',
             kv:  k ? k.textContent : '',
             txt: (c.textContent || '').toLowerCase() };
  });
  var cur = -1;

  /* ---------- lightbox ---------- */
  var lb = document.getElementById('lb');
  var img = document.getElementById('lb-img');
  var noEl = document.getElementById('lb-no');
  var ttlEl = document.getElementById('lb-ttl');
  var kvEl = document.getElementById('lb-kv');
  var capEl = document.getElementById('lb-cap');
  var dl = document.getElementById('lb-dl');

  function show(i){
    if (i < 0) i = cards.length - 1;
    if (i >= cards.length) i = 0;
    cur = i;
    var m = meta[i];
    img.src = FULL[m.idx];
    noEl.textContent = m.no;
    ttlEl.textContent = m.ttl;
    kvEl.textContent = m.kv;
    var cap = cards[i].querySelector('.cap');
    capEl.innerHTML = cap ? cap.innerHTML : '';
    dl.download = m.no + '_' + (m.ttl.replace(/[\\/:*?"<>|\s]+/g,'_')) + '.png';
    dl.href = FULL[m.idx];
    lb.classList.add('on');
    document.body.style.overflow = 'hidden';
    var st = lb.querySelector('.stage'); st.scrollTop = 0; st.scrollLeft = 0;
  }
  function hide(){
    lb.classList.remove('on');
    img.src = '';
    document.body.style.overflow = '';
  }
  cards.forEach(function(c, i){
    c.addEventListener('click', function(){ show(i); });
  });
  document.getElementById('lb-close').onclick = hide;
  document.getElementById('lb-prev').onclick = function(e){ e.stopPropagation(); show(cur-1); };
  document.getElementById('lb-next').onclick = function(e){ e.stopPropagation(); show(cur+1); };
  lb.addEventListener('click', function(e){ if (e.target === lb || e.target.className === 'stage') hide(); });
  document.addEventListener('keydown', function(e){
    if (!lb.classList.contains('on')) return;
    if (e.key === 'Escape') hide();
    else if (e.key === 'ArrowLeft')  show(cur-1);
    else if (e.key === 'ArrowRight') show(cur+1);
  });

  /* ---------- filter ---------- */
  var box = document.getElementById('q');
  var cnt = document.getElementById('qcount');
  box.addEventListener('input', function(){
    var q = box.value.trim().toLowerCase();
    var n = 0;
    cards.forEach(function(c, i){
      var hit = !q || meta[i].txt.indexOf(q) >= 0;
      c.classList.toggle('hidden', !hit);
      if (hit) n++;
    });
    [].slice.call(document.querySelectorAll('.grp-title')).forEach(function(g){
      var any = false, el = g.nextElementSibling;
      while (el && el.classList && el.classList.contains('grid')) {
        if ([].slice.call(el.children).some(function(x){ return !x.classList.contains('hidden'); })) any = true;
        el = el.nextElementSibling;
      }
      g.classList.toggle('hidden', !any);
    });
    cnt.textContent = q ? (n + ' / ' + cards.length + ' 张匹配') : (cards.length + ' 张');
  });
  cnt.textContent = cards.length + ' 张';

  /* ---------- expand / collapse all ---------- */
  var exp = false;
  document.getElementById('expand').onclick = function(){
    exp = !exp;
    document.querySelectorAll('.thumb img').forEach(function(im){
      im.style.maxHeight = exp ? 'none' : '';
    });
    this.textContent = exp ? '收起缩略图' : '展开缩略图';
  };
})();
"""

TOOLBAR = """
<div class="toolbar">
  <input id="q" type="search" placeholder="筛选：图号 / 标题 / 关键词（如 RMSD、H 键、Mn、DCCM）…" autocomplete="off">
  <span class="count" id="qcount"></span>
  <button id="expand" class="toc-grp" style="border:1px solid #2A3648;cursor:pointer;color:var(--ink)">展开缩略图</button>
  <span class="hint">点击任意图片 → 原图无损放大  ·  键盘 ←/→ 翻页  ·  Esc 关闭</span>
</div>
"""

LIGHTBOX = """
<div id="lb">
  <div class="bar">
    <span class="no" id="lb-no"></span>
    <span class="ttl" id="lb-ttl"></span>
    <span class="sp"></span>
    <button id="lb-prev">← 上一张</button>
    <button id="lb-next">下一张 →</button>
    <a id="lb-dl" download><button>下载原图</button></a>
    <button id="lb-close">关闭 ✕</button>
  </div>
  <div class="stage"><img id="lb-img" alt=""></div>
  <div class="foot"><span id="lb-kv"></span><br><span id="lb-cap"></span></div>
</div>
"""


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=None, help="figures_pretty 目录（含 index.html + png）")
    ap.add_argument("--out", default=None, help="输出单文件 HTML 路径")
    ap.add_argument("--thumb-width", type=int, default=800)
    ap.add_argument("--thumb-quality", type=int, default=82)
    ap.add_argument("--full", choices=["png", "jpeg"], default="png")
    ap.add_argument("--full-quality", type=int, default=92)
    a = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    fdir = a.dir or os.path.join(here, "results", "WT__S1", "figures_pretty")
    fdir = os.path.abspath(fdir)
    idx_html = os.path.join(fdir, "index.html")
    if not os.path.exists(idx_html):
        sys.exit("[err] 找不到 " + idx_html)

    html = open(idx_html, encoding="utf-8").read()

    # 1) 收集卡片顺序 (按 id=figNN 出现顺序)
    cards = re.findall(r'<a class="card" id="(fig\d+)" href="([^"]+\.png)"', html)
    if not cards:
        sys.exit("[err] index.html 中未解析到卡片，请检查结构")
    print(f"[1/4] 解析 index.html: {len(cards)} 张图")

    # 2) 缩略图 + 原图 base64
    thumbs, fulls = {}, {}
    for i, (fid, png) in enumerate(cards):
        p = os.path.join(fdir, png)
        if not os.path.exists(p):
            sys.exit(f"[err] 缺图: {p}")
        thumbs[png] = b64_thumb(p, a.thumb_width, a.thumb_quality)
        fulls[fid] = b64_full(p, a.full, a.full_quality)
        print(f"      [{i+1:2d}/{len(cards)}] {png}"
              f"  thumb {len(thumbs[png])/1024:6.0f} KB"
              f"  full {len(fulls[fid])/1024:7.0f} KB")
    print("[2/4] base64 编码完成")

    # 3) 改写 DOM
    #    a. 卡片 a -> div (JS 驱动灯箱)
    for fid, png in cards:
        html = html.replace(
            f'<a class="card" id="{fid}" href="{png}" target="_blank">',
            f'<div class="card" id="{fid}" data-idx="{cards.index((fid, png))}">', 1)
    #    b. 图片 src -> 缩略图 data URI
    def _sub(m):
        png = m.group(1)
        return f'src="{thumbs[png]}"'
    html, n_img = re.subn(r'src="([^"]+\.png)"', _sub, html)
    #    c. 关闭 </a> —— <main> 内剩余的 </a> 全是卡片闭合标签（header/toc 的 <a> 在 main 之外）
    i0, i1 = html.index("<main>"), html.index("</main>")
    main_html = html[i0:i1]
    n_cls = main_html.count("</a>")
    html = html[:i0] + main_html.replace("</a>", "</div>") + html[i1:]
    leftover = len(re.findall(r'<div class="card"', html)) - html.count('<div class="card"')
    print(f"[3/4] DOM 改写: {n_img} 张 img 换缩略图, {n_cls} 个卡片闭合标签修复")
    if main_html.count('<div class="card"') != n_cls:
        print(f"      [warn] 卡片开标签 {main_html.count(chr(60)+'div class=' + chr(34) + 'card' + chr(34))} != 闭标签 {n_cls}，请人工核对")

    # 4) 注入 CSS / 工具栏 / 灯箱 / JS
    html = html.replace("</style>", EXTRA_CSS + "\n</style>", 1)
    html = html.replace('<div class="toc">', TOOLBAR + '<div class="toc">', 1)
    html = html.replace("</body>", LIGHTBOX + "\n<script>\n"
                        + "window.__FULL__ = [\n"
                        + ",\n".join('"%s"' % fulls[fid] for fid, _ in cards)
                        + "\n];\n</script>\n<script>\n" + EXTRA_JS + "\n</script>\n</body>", 1)

    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    html = html.replace("<header>",
        f'<!-- 自包含单文件图集 · 生成于 {stamp} · '
        f'{len(cards)} 张图 · 无外部依赖，可直接分享 -->\n<header>', 1)
    html = html.replace('</h1>',
        '</h1>\n  <div class="meta" style="margin-top:6px">'
        '<span style="color:#E8623C">单文件版</span>：30 张图全部内嵌 base64，'
        '无需文件夹、断网可用 · 缩略图 = 轻量 JPEG，点击查看/下载原图 PNG'
        f' · 生成于 {stamp}</div>', 1)

    out = a.out or os.path.join(os.path.dirname(fdir),
                                os.path.basename(os.path.dirname(fdir)) + "_figures_standalone.html")
    out = os.path.abspath(out)
    with open(out, "w", encoding="utf-8") as f:
        f.write(html)

    size_mb = os.path.getsize(out) / 1024 / 1024
    print(f"[4/4] 输出: {out}")
    print(f"      文件大小 {size_mb:.1f} MB · {len(cards)} 张图 · 完全自包含")
    print("[done] 双击即可打开；可直接发送给他人，无需打包文件夹。")


if __name__ == "__main__":
    main()
