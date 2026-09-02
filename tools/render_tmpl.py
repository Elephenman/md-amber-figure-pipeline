#!/usr/bin/env python3
"""render_tmpl.py — 用 python 做占位符替换（避免 sed 对 mask 中 | & 特殊字符的兼容问题）
用法: render_tmpl.py <模板> <输出> KEY=VAL [KEY=VAL ...]
模板中 {KEY} 占位符会被替换为 VAL。"""
import sys, os

def main():
    if len(sys.argv) < 3:
        sys.stderr.write("用法: render_tmpl.py <模板> <输出> KEY=VAL [...]\n")
        return 2
    tpl, out = sys.argv[1], sys.argv[2]
    with open(tpl, 'r') as f:
        text = f.read()
    for kv in sys.argv[3:]:
        k, _, v = kv.partition('=')
        text = text.replace('{' + k + '}', v)
    d = os.path.dirname(out)
    if d:  # 相对路径 dirname 为空，跳过 makedirs
        os.makedirs(d, exist_ok=True)
    with open(out, 'w') as f:
        f.write(text)
    return 0

if __name__ == '__main__':
    sys.exit(main())
