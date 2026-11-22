import json, sys, glob, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.ticker
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.abspath(__file__)))
SITE, GAUSS = "#2a78d6", "#eb6834"
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "font.size": 10,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "lines.linewidth": 2, "lines.markersize": 7})
cases = [("uniform", 32), ("uniform", 64), ("uniform", 128), ("wilson", 32)]
titles = {("uniform", 32): "uniform bath, 32 orbitals", ("uniform", 64): "uniform bath, 64 orbitals",
          ("uniform", 128): "uniform bath, 128 orbitals", ("wilson", 32): "Wilson chain, 32 orbitals"}
data = {}
for c in cases:
    f = glob.glob(os.path.join("data", "cmp_%s_n%d_*.json" % c))
    if f: data[c] = json.load(open(f[0]))
fig, ax = plt.subplots(2, len(data), figsize=(4.2*len(data), 7.4), squeeze=False)
for k, (c, d) in enumerate(data.items()):
    runs = d["runs"]; eref = min(r["E"] for r in runs)
    for basis, col, mk, lab in [("site", SITE, "o", "site basis"), ("gaussian", GAUSS, "s", "Gaussian basis")]:
        rr = [r for r in runs if r["basis"] == basis]
        m = np.array([r["maxm"] for r in rr]); err = np.array([r["E"]-eref for r in rr]); t = np.array([r["t"] for r in rr])
        ok = err > 0
        ax[0, k].loglog(m[ok], err[ok], mk+"-", color=col, label=lab)
        ax[1, k].loglog(t[ok], err[ok], mk+"-", color=col, label=lab)
    ax[0, k].set_title(titles[c], color=INK, fontsize=11)
    ax[0, k].set_xlabel("bond dimension maxm"); ax[1, k].set_xlabel("wall time of the solve (s)")
    ax[0, k].set_xticks([8, 16, 32, 64, 128]); ax[0, k].set_xticklabels(["8", "16", "32", "64", "128"])
    ax[0, k].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax[1, k].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax[0, 0].set_ylabel("E - E_best (units of t)"); ax[1, 0].set_ylabel("E - E_best (units of t)")
ax[0, 0].legend(frameon=False, labelcolor=INK)
fig.suptitle("Anderson impurity, U=4 (uniform bath) and U=0.5 (Wilson chain): DMRG energy error, v3, 12 sweeps, cold start; time includes the MPO build", color=INK, fontsize=12)
fig.tight_layout(); fig.savefig("compare_energy.png", dpi=130)

fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
c = ("uniform", 128) if ("uniform", 128) in data else ("uniform", 64)
for basis, col, lab in [("site", SITE, "site basis"), ("gaussian", GAUSS, "Gaussian basis")]:
    r = max([r for r in data[c]["runs"] if r["basis"] == basis], key=lambda r: r["maxm"])
    S = np.array(r["S"]); x = np.arange(1, len(S)+1)/2
    ax[0].plot(x, S, color=col, linewidth=1.5, label="%s, maxm=%d" % (lab, r["maxm"]))
ax[0].set_xlabel("cut position (orbitals from the impurity)"); ax[0].set_ylabel("entanglement entropy")
ax[0].set_title(titles[c], color=INK, fontsize=11); ax[0].legend(frameon=False, labelcolor=INK)
ns = sorted(n for (case, n) in data if case == "uniform")
for basis, col, mk, lab in [("site", SITE, "o", "site basis"), ("gaussian", GAUSS, "s", "Gaussian basis")]:
    ax[1].plot(ns, [data[("uniform", n)]["runs"][[r["basis"] for r in data[("uniform", n)]["runs"]].index(basis)]["mpo"] for n in ns], mk+"-", color=col, label=lab)
ax[1].set_xscale("log", base=2); ax[1].set_xticks(ns); ax[1].set_xticklabels([str(n) for n in ns])
ax[1].set_xlabel("orbitals (impurity + bath)"); ax[1].set_ylabel("MPO bond dimension of H")
ax[1].set_title("cost of the rotation: MPO of H", color=INK, fontsize=11); ax[1].legend(frameon=False, labelcolor=INK)
fig.tight_layout(); fig.savefig("compare_entropy_mpo.png", dpi=130)
for c, d in data.items():
    eref = min(r["E"] for r in d["runs"])
    print(c, "E_best=%.10f" % eref)
    for r in d["runs"]:
        print("  %-8s maxm=%4d err=%.2e t=%7.1fs mpo=%d Smax=%.3f docc=%.6f" % (r["basis"], r["maxm"], r["E"]-eref, r["t"], r["mpo"], max(r["S"]), r["docc"]))
