"""Timings and profiles of the Kondo calls of docs/kondo_speedup_plan.md.

One case is one public call. For every code tree given with --src, and for
every repeat, the call runs in a fresh Python process (the order of the
trees alternates between repeats), so two versions of the code are timed on
the same node in the same job. One more run per tree is under cProfile and
gives the share of the functions the plan names. The result goes to one
JSON file: the CPU model, the wall times, the peak resident memory, the
profile shares and the largest difference of the curve between trees.

    python run_case.py --list
    python run_case.py --case ed_T0_n10 --src /path/to/src --json out.json
    python run_case.py --case ed_T0_n10 --src old/src --src new/src --reps 3

Pin BLAS threads before running (MKL_NUM_THREADS=1 OMP_NUM_THREADS=1).
"""
import argparse
import json
import os
import platform
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

CHILD = r'''
import cProfile, io, json, pstats, resource, sys, time, warnings
import numpy as np
warnings.simplefilter("ignore")
sys.path.insert(0, sys.argv[1])
sys.path.insert(0, sys.argv[2])
import cases
case = cases.CASES[sys.argv[3]]
profile = sys.argv[4] == "1"
setup = case["setup"]()
prof = cProfile.Profile() if profile else None
t0 = time.perf_counter()
if prof: prof.enable()
out = case["call"](setup)
if prof: prof.disable()
t = time.perf_counter() - t0
res = dict(t=t, peak_kb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
           curve=np.real(np.asarray(out, dtype=complex)).ravel().tolist())
if prof:
    st = pstats.Stats(prof)
    share = {}
    for (fname, line, name), (cc, nc, tt, ct, callers) in st.stats.items():
        for key in case["watch"]:
            if key == name or key == "%s:%s" % (fname.split("/")[-1], name):
                share[key] = max(share.get(key, 0.), ct)
    res["profile_s"] = share
    s = io.StringIO()
    pstats.Stats(prof, stream=s).sort_stats("cumulative").print_stats(25)
    res["top"] = s.getvalue().splitlines()[:60]
print("RESULT " + json.dumps(res))
'''


def cpu_model():
    try:
        for line in open("/proc/cpuinfo"):
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor()


def run(src, case, profile):
    env = dict(os.environ, PYTHONPATH=src)
    out = subprocess.run([sys.executable, "-s", "-c", CHILD, src, HERE, case,
                          "1" if profile else "0"],
                         capture_output=True, text=True, env=env)
    for line in out.stdout.splitlines():
        if line.startswith("RESULT "):
            return json.loads(line[7:])
    raise RuntimeError("case %s on %s failed:\n%s" % (case, src, out.stderr[-3000:]))


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--case")
    p.add_argument("--src", action="append", default=[])
    p.add_argument("--reps", type=int, default=None)
    p.add_argument("--no-profile", action="store_true")
    p.add_argument("--json")
    p.add_argument("--list", action="store_true")
    p.add_argument("--inline", action="store_true",
                   help="run the case once in this process, for a smoke test "
                        "(one Python process; not for timing two trees)")
    a = p.parse_args()
    # cases.py imports dmrgpy: from the first tree given, or from this
    # checkout, never from whatever an installed copy would give
    src0 = os.path.abspath(a.src[0]) if a.src else os.path.join(HERE, "..", "..", "src")
    sys.path.insert(0, src0)
    sys.path.insert(0, HERE)
    import cases
    if a.list:
        for i, (k, c) in enumerate(cases.CASES.items()):
            print(i, k, "--", c["what"])
        return
    if a.case.isdigit():
        a.case = list(cases.CASES)[int(a.case)]
    case = cases.CASES[a.case]
    if a.inline:
        import time
        import resource
        setup = case["setup"]()
        t0 = time.perf_counter()
        case["call"](setup)
        print("%s inline: %.2f s, peak %.2f GB" % (a.case, time.perf_counter() - t0,
              resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6))
        return
    reps = case["reps"] if a.reps is None else a.reps
    srcs = [os.path.abspath(s) for s in a.src]
    report = dict(case=a.case, what=case["what"], cpu=cpu_model(),
                  host=platform.node(), srcs=srcs, reps=reps, runs={})
    for s in srcs:
        report["runs"][s] = dict(t=[], peak_kb=[])
    curves = {}
    for r in range(reps):
        for s in (srcs if r % 2 == 0 else srcs[::-1]):
            res = run(s, a.case, False)
            report["runs"][s]["t"].append(res["t"])
            report["runs"][s]["peak_kb"].append(res["peak_kb"])
            curves.setdefault(s, res["curve"])
            print("%s rep %d %s: %.2f s, %.2f GB" % (a.case, r, s, res["t"],
                  res["peak_kb"]/1e6), flush=True)
    if not a.no_profile:
        for s in srcs:
            res = run(s, a.case, True)
            report["runs"][s]["profile_total_s"] = res["t"]
            report["runs"][s]["profile_s"] = res["profile_s"]
            report["runs"][s]["top"] = res["top"]
            curves.setdefault(s, res["curve"])
            print("%s profile %s: %.2f s, %s" % (a.case, s, res["t"],
                  {k: round(v, 2) for k, v in res["profile_s"].items()}), flush=True)
    if len(srcs) > 1:
        import numpy as np
        ref = np.asarray(curves[srcs[0]])
        scale = max(np.max(np.abs(ref)), 1e-300)
        report["max_rel_diff"] = {s: float(np.max(np.abs(np.asarray(curves[s]) - ref))/scale)
                                  for s in srcs[1:]}
    if a.json:
        with open(a.json, "w") as f:
            json.dump(report, f, indent=1)


if __name__ == "__main__":
    main()
