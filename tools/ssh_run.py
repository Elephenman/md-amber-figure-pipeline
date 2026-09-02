#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ssh_run.py — CHPC 通用远程执行 + 文件拉取/推送
用法:
  CHPC_PASS=xxx python ssh_run.py "cmd"          # 执行远程命令
  CHPC_PASS=xxx python ssh_run.py pull a b c     # 从远端拉文件(自动下载到 cwd 下 mirror/)
  CHPC_PASS=xxx python ssh_run.py push l1 r1     # 上传本地文件到远端绝对路径
"""
import os
import sys
import posixpath
import paramiko

HOST, PORT, USER = "10.202.94.52", 20009, "u22607007"


def connect():
    pw = os.environ.get("CHPC_PASS")
    if not pw:
        raise SystemExit("ERROR: 请设置 CHPC_PASS")
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(HOST, port=PORT, username=USER, password=pw, timeout=30,
              look_for_keys=False, allow_agent=False)
    return c


def run(c, cmd, timeout=900):
    stdin, stdout, stderr = c.exec_command(cmd, timeout=timeout)
    out = stdout.read().decode(errors="replace")
    err = stderr.read().decode(errors="replace")
    rc = stdout.channel.recv_exit_status()
    if out.strip():
        print(out)
    if err.strip():
        print("[stderr]", err[:3000])
    print(f"==== exit {rc} ====")
    return rc


def pull(c, remote_files):
    sftp = c.open_sftp()
    for rp in remote_files:
        fn = posixpath.basename(rp)
        local = os.path.join("mirror", fn)
        os.makedirs("mirror", exist_ok=True)
        sftp.get(rp, local)
        print(f"[pull] {rp} -> {local}")
    sftp.close()


def push(c, local, remote):
    sftp = c.open_sftp()
    rdir = posixpath.dirname(remote)
    acc = "/" if rdir.startswith("/") else ""
    for seg in rdir.split("/"):
        if not seg:
            continue
        acc = posixpath.join(acc, seg)
        try:
            sftp.stat(acc)
        except IOError:
            sftp.mkdir(acc)
    sftp.put(local, remote)
    print(f"[push] {local} -> {remote}")
    sftp.close()


def main():
    args = sys.argv[1:]
    if not args:
        raise SystemExit(__doc__)
    c = connect()
    try:
        if args[0] == "pull":
            pull(c, args[1:])
        elif args[0] == "push":
            push(c, args[1], args[2])
        else:
            run(c, " ".join(args))
    finally:
        c.close()


if __name__ == "__main__":
    main()
