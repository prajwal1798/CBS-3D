#!/usr/bin/env python3
# ================================
# FILE      : scratch_inventory.py
# PURPOSE   : READ-ONLY MAP OF EVERYTHING CBS3D-RELATED IN YOUR SUNBIRD SCRATCH
#             AND HOME: CODE COPIES, EXECUTABLES, CASES, PARTITION SETS, RUNS,
#             VTU OUTPUT, RESTART CHECKPOINTS, JOB SCRIPTS, LOGS, JOB HISTORY.
#             WRITES ONE TEXT FILE TO ATTACH IN A CHAT.
#
# USAGE     : python3 scratch_inventory.py              (on the Sunbird login node)
#             python3 scratch_inventory.py --root /scratch/s.2337862 --home ~
#
# SAFETY    : NEVER WRITES, MOVES OR DELETES ANYTHING EXCEPT THE REPORT FILE.
# STATUS    : PYTHON 3.6 COMPATIBLE. A FULL WALK OF A LARGE SCRATCH TAKES MINUTES.
# ================================

from __future__ import print_function

import argparse
import datetime
import os
import re
import subprocess
import sys
import time

STEP = re.compile(r"_step_(\d+)")
RANK_DIR = re.compile(r"^rank_\d{4}$")
CODE_EXE = ("cbs3d_parallel", "cbs3d_serial", "cbs3dpp_si", "cbs3dpp_si_petsc")
SCRIPT_EXT = (".slurm", ".sh", ".sbatch")
KEY_EXT = (".geo", ".msh")


# ================================
# SMALL HELPERS
# ================================
def run(cmd, cwd=None, timeout=120):
    try:
        p = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             universal_newlines=True)
        out, err = p.communicate(timeout=timeout)
        return p.returncode, out.strip()
    except Exception as e:  # noqa
        return 1, str(e)


def human(n):
    for unit in ("B", "K", "M", "G", "T"):
        if n < 1024.0:
            return "%.1f%s" % (n, unit)
        n /= 1024.0
    return "%.1fP" % n


def stamp(t):
    return datetime.datetime.fromtimestamp(t).strftime("%Y-%m-%d %H:%M")


def head_lines(path, n):
    out = []
    try:
        with open(path, "r", errors="replace") as f:
            for i, line in enumerate(f):
                if i >= n:
                    break
                out.append(line.rstrip("\n"))
    except Exception:
        pass
    return out


def tail_lines(path, n, max_bytes=200000):
    try:
        size = os.path.getsize(path)
        with open(path, "rb") as f:
            f.seek(max(0, size - max_bytes))
            data = f.read().decode("utf-8", "replace")
        return data.splitlines()[-n:]
    except Exception:
        return []


def par_values(path):
    # POSITIONAL .par: LABEL LINE, DATA LINE. RETURN THE BLOCKS WE CARE ABOUT.
    lines = head_lines(path, 80)
    blocks = {}
    for i in range(0, len(lines) - 1, 2):
        label = lines[i].strip().lower()
        data = lines[i + 1].strip()
        if i == 0:
            blocks["solver_opt"] = data
        elif "restart" in label:
            blocks["restart"] = data
        elif "energy" in label:
            blocks["energy"] = data
        elif label.startswith("ntime"):
            blocks["ntime transient dtfixed dtfix iwrite"] = data
        elif "output/monitor" in label or label.startswith("residual_log_enabled"):
            blocks["output/monitor (line 44)"] = data
        elif "turbulence_on" in label:
            blocks["turbulence"] = data
    return blocks


# ================================
# ONE WALK OVER THE TREE, COLLECTING EVERYTHING
# ================================
def walk(root, skip_hidden=True):
    info = {
        "dirs": {},          # path -> [n_files, bytes, newest_mtime]
        "git": [],           # repo roots
        "exe": [],           # (path, size, mtime)
        "serial_cases": [],  # .plt not inside rank dirs
        "partition_sets": [],
        "outputs": {},       # dir -> {steps:set, pieces, pvtu, csv:[...], restart:[]}
        "manifests": [],
        "scripts": [],
        "logs": [],
        "meshes": [],
        "pars": [],
    }
    for dirpath, dirnames, filenames in os.walk(root):
        if skip_hidden:
            dirnames[:] = [d for d in dirnames if not d.startswith(".") or d == ".git"]
        if ".git" in dirnames:
            info["git"].append(dirpath)
            dirnames.remove(".git")
        if "CMakeFiles" in dirnames:
            dirnames.remove("CMakeFiles")
        if any(RANK_DIR.match(d) for d in dirnames):
            ranks = sorted(d for d in dirnames if RANK_DIR.match(d))
            info["partition_sets"].append((dirpath, ranks))
        nbytes = 0
        newest = 0.0
        out = None
        for name in filenames:
            path = os.path.join(dirpath, name)
            try:
                st = os.lstat(path)
            except OSError:
                continue
            nbytes += st.st_size
            newest = max(newest, st.st_mtime)
            low = name.lower()
            if name in CODE_EXE and os.access(path, os.X_OK):
                info["exe"].append((path, st.st_size, st.st_mtime))
            elif low.endswith(".plt") and "_rank_" not in name:
                info["serial_cases"].append(path)
            elif low.endswith(KEY_EXT):
                info["meshes"].append((path, st.st_size, st.st_mtime))
            elif low.endswith(".par") and "_rank_" not in name:
                info["pars"].append(path)
            elif name == "restart_manifest.txt":
                info["manifests"].append(path)
            elif low.endswith(SCRIPT_EXT):
                info["scripts"].append((path, st.st_mtime))
            elif low.endswith((".out", ".err", ".log")) and st.st_size > 0:
                info["logs"].append((path, st.st_size, st.st_mtime))
            if low.endswith((".vtu", ".pvtu", ".csv", ".pvd", ".cbsrst")):
                if out is None:
                    out = info["outputs"].setdefault(dirpath, {"steps": set(), "pieces": 0, "pvtu": 0,
                                                             "csv": [], "bytes": 0, "rst": 0})
                out["bytes"] += st.st_size
                m = STEP.search(name)
                if m and low.endswith((".vtu", ".pvtu", ".cbsrst")):
                    out["steps"].add(int(m.group(1)))
                if low.endswith(".vtu"):
                    out["pieces"] += 1
                elif low.endswith(".pvtu"):
                    out["pvtu"] += 1
                elif low.endswith(".csv"):
                    out["csv"].append(path)
                elif low.endswith(".cbsrst"):
                    out["rst"] += 1
        info["dirs"][dirpath] = [len(filenames), nbytes, newest]
    return info


def tree_sizes(info, root):
    # ACCUMULATE BYTES AND FILES UP THE TREE
    total = {}
    for d, (nf, nb, nt) in info["dirs"].items():
        p = d
        while True:
            t = total.setdefault(p, [0, 0, 0.0])
            t[0] += nf
            t[1] += nb
            t[2] = max(t[2], nt)
            if p == root or len(p) <= len(root):
                break
            p = os.path.dirname(p)
    return total


# ================================
# REPORT SECTIONS
# ================================
def section(out, title):
    out.append("")
    out.append("=" * 100)
    out.append(title)
    out.append("=" * 100)


def report_tree(out, root, total, max_depth):
    section(out, "TREE OF %s (depth <= %d; rank_XXXX folders collapsed)" % (root, max_depth))
    out.append("%10s %9s  %-16s  %s" % ("SIZE", "FILES", "NEWEST", "FOLDER"))
    base = root.rstrip("/").count("/")
    for d in sorted(total):
        depth = d.rstrip("/").count("/") - base
        if depth > max_depth:
            continue
        name = os.path.basename(d)
        if RANK_DIR.match(name) and name != "rank_0000":
            continue
        nf, nb, nt = total[d]
        label = d if depth == 0 else "  " * depth + name + "/"
        if name == "rank_0000":
            label += "   (+ other rank_XXXX)"
        out.append("%10s %9d  %-16s  %s" % (human(nb), nf, stamp(nt) if nt else "-", label))


def report_code(out, info):
    section(out, "CODE COPIES (git repositories) AND EXECUTABLES")
    for repo in sorted(info["git"]):
        _, branch = run(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=repo)
        _, last = run(["git", "log", "-1", "--format=%h %ad %s", "--date=short"], cwd=repo)
        _, dirty = run(["git", "status", "--porcelain"], cwd=repo)
        _, remote = run(["git", "remote", "-v"], cwd=repo)
        out.append("REPO   %s" % repo)
        out.append("  branch : %s" % branch)
        out.append("  commit : %s" % last)
        out.append("  local changes : %d files" % (len(dirty.splitlines()) if dirty else 0))
        out.append("  remote : %s" % (remote.splitlines()[0] if remote else "none"))
    out.append("")
    for path, size, mt in sorted(info["exe"], key=lambda x: -x[2]):
        out.append("EXE    %-16s %8s  %s" % (stamp(mt), human(size), path))


def report_cases(out, info):
    section(out, "SERIAL CASES (.plt not in rank folders) WITH HEADER")
    for p in sorted(info["serial_cases"]):
        h = head_lines(p, 1)
        out.append("%-90s  nelem npoin nboun = %s" % (p, h[0].strip() if h else "?"))
    section(out, "PARTITION SETS (folders holding rank_XXXX)")
    for d, ranks in sorted(info["partition_sets"]):
        first = os.path.join(d, ranks[0])
        files = sorted(os.listdir(first)) if os.path.isdir(first) else []
        case = "?"
        for f in files:
            m = re.match(r"(.+)_rank_\d{4}\.plt$", f)
            if m:
                case = m.group(1)
        mpi = [l for l in head_lines(os.path.join(first, "%s_rank_0000.mpi" % case), 4) if l.startswith("GLOBAL")]
        out.append("SET    %s" % d)
        out.append("  ranks=%d  case=%s  %s" % (len(ranks), case, mpi[0] if mpi else ""))
        par = os.path.join(first, "%s_rank_0000.par" % case)
        for k, v in par_values(par).items():
            out.append("  par %-40s : %s" % (k, v))
    section(out, ".par FILES OUTSIDE RANK FOLDERS (KEY BLOCKS)")
    for p in sorted(info["pars"]):
        out.append(p)
        for k, v in par_values(p).items():
            out.append("    %-40s : %s" % (k, v))


def report_outputs(out, info):
    section(out, "OUTPUT FOLDERS (VTU / PVTU / CSV / NATIVE RESTART)")
    for d in sorted(info["outputs"]):
        o = info["outputs"][d]
        steps = sorted(o["steps"])
        out.append("OUT    %s" % d)
        out.append("  size %s | pieces(.vtu) %d | .pvtu %d | .cbsrst %d" % (human(o["bytes"]), o["pieces"], o["pvtu"], o["rst"]))
        if steps:
            out.append("  steps: %d distinct, first %d, last %d" % (len(steps), steps[0], steps[-1]))
        for c in o["csv"]:
            last = tail_lines(c, 1)
            hdr = head_lines(c, 1)
            out.append("  csv %s" % os.path.basename(c))
            if hdr:
                out.append("     header: %s" % hdr[0][:300])
            if last:
                out.append("     last  : %s" % last[0][:300])
    section(out, "NATIVE RESTART MANIFESTS")
    for m in sorted(info["manifests"]):
        out.append("MANIFEST %s" % m)
        for line in head_lines(m, 20):
            out.append("    " + line)


def report_scripts(out, info):
    section(out, "JOB SCRIPTS (.slurm .sh) - key lines")
    keys = re.compile(r"(#SBATCH --(job-name|ntasks|nodes|time|cpus-per-task)|CBS3D_[A-Z_]+=|EXE=|CASE=|ROOT=|PARTITION_ROOT=|SOURCE_PART=|srun|mpirun|sbatch )")
    for path, mt in sorted(info["scripts"], key=lambda x: -x[1]):
        out.append("SCRIPT %s  (%s)" % (path, stamp(mt)))
        for line in head_lines(path, 400):
            s = line.strip()
            if keys.search(s) and not s.startswith("echo") and len(s) < 220:
                out.append("    " + s)


def report_logs(out, info, n_logs):
    section(out, "NEWEST %d LOG FILES - last informative lines" % n_logs)
    pat = re.compile(r"(stop reason|Run complete|Restart loaded|Continuation|Iteration \d+|it=\d+|FATAL|ERROR|CANCELLED|TIME LIMIT|RelU|Stop reason|Last iteration)", re.I)
    for path, size, mt in sorted(info["logs"], key=lambda x: -x[2])[:n_logs]:
        out.append("LOG    %s  (%s, %s)" % (path, stamp(mt), human(size)))
        hits = [l for l in tail_lines(path, 400) if pat.search(l)]
        for l in hits[-6:]:
            out.append("    " + re.sub(r"\x1b\[[0-9;]*m", "", l)[:220])


def report_meshes(out, info):
    section(out, ".geo AND .msh FILES")
    for path, size, mt in sorted(info["meshes"]):
        out.append("%-16s %9s  %s" % (stamp(mt), human(size), path))


def report_sacct(out, days):
    section(out, "SLURM JOB HISTORY (last %d days)" % days)
    start = (datetime.date.today() - datetime.timedelta(days=days)).isoformat()
    code, text = run(["sacct", "-X", "-S", start, "-o",
                      "JobID,JobName%28,State,Elapsed,AllocCPUS,End,WorkDir%110"], timeout=120)
    out.extend(text.splitlines() if text else ["(sacct not available: code %d)" % code])


# ================================
# MAIN
# ================================
def main():
    user = os.environ.get("USER", "")
    ap = argparse.ArgumentParser(description="Read-only CBS3D inventory of Sunbird scratch and home.")
    ap.add_argument("--root", default="/scratch/%s" % user)
    ap.add_argument("--home", default=os.path.expanduser("~"))
    ap.add_argument("--depth", type=int, default=4, help="tree depth to print")
    ap.add_argument("--logs", type=int, default=40, help="newest log files to summarise")
    ap.add_argument("--days", type=int, default=150, help="days of sacct history")
    ap.add_argument("--out", default=os.path.expanduser("~/scratch_inventory_%s.txt" % datetime.date.today().isoformat()))
    a = ap.parse_args()

    t0 = time.time()
    out = ["CBS3D SCRATCH INVENTORY", "date   : %s" % datetime.datetime.now().isoformat(),
           "host   : %s" % os.uname()[1], "scratch: %s" % a.root, "home   : %s" % a.home]

    print("walking %s ..." % a.root)
    info = walk(a.root)
    print("walking %s ..." % a.home)
    hinfo = walk(a.home)
    for k in ("git", "exe", "serial_cases", "partition_sets", "manifests", "scripts", "logs", "meshes", "pars"):
        info[k].extend(hinfo[k])
    info["outputs"].update(hinfo["outputs"])

    report_tree(out, a.root, tree_sizes(info, a.root), a.depth)
    report_tree(out, a.home, tree_sizes(hinfo, a.home), 2)
    report_code(out, info)
    report_cases(out, info)
    report_outputs(out, info)
    report_meshes(out, info)
    report_scripts(out, info)
    report_logs(out, info, a.logs)
    report_sacct(out, a.days)
    out.append("")
    out.append("finished in %.0f s" % (time.time() - t0))

    with open(a.out, "w") as f:
        f.write("\n".join(out) + "\n")
    print("REPORT WRITTEN: %s (%s)" % (a.out, human(os.path.getsize(a.out))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
