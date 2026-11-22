"""Energy error against the best run, per basis and maxm, and the MPO bond
dimension, from the data/ref_*.json files of reference.py."""
import glob, json, os
here = os.path.dirname(os.path.abspath(__file__))
for f in sorted(glob.glob(os.path.join(here, "data", "ref_*.json"))):
    d = json.load(open(f))
    best = min(r["E"] for r in d["runs"])
    print("\n%s n=%d U=%g ed=%g nsweeps=%d  levels %s  E_best=%.10f" % (d["case"], d["n"], d["U"], d["ed"], d.get("nsweeps", 12),
          {k: round(v, 4) for k, v in d.get("levels", {}).items()}, best))
    if "spin_mixing" in d: print("  spin mixing dropped %.1e" % d["spin_mixing"])
    rows = {}
    for r in d["runs"]:
        rows.setdefault(r["basis"], {})[r["maxm"]] = (r["E"]-best, r["mpo"])
    ms = sorted({r["maxm"] for r in d["runs"]})
    print("  %-14s mpo " % "basis" + " ".join("%9d" % m for m in ms))
    for b, row in rows.items():
        mpo = list(row.values())[0][1]
        print("  %-14s %3d " % (b, mpo) + " ".join("%9.1e" % row[m][0] if m in row else "%9s" % "-" for m in ms))
