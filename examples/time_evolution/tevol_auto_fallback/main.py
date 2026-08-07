# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Demonstrates tevol_method="AUTO" (timedependent._tebd_or_tdvp()): a
# script that doesn't know in advance whether its own Hamiltonian is
# strictly nearest-neighbor -- e.g. because it sweeps a hopping range, or
# is reused across several models -- can just set tevol_method="AUTO"
# once and let it try the cheaper "TEBD" propagator first, transparently
# retrying with "TDVP" only when the Hamiltonian actually needs it,
# instead of the caller hand-picking the method per case.
#
# Base model: nearest-neighbor hopping plus a staggered onsite field (the
# field breaks the would-be degeneracy of a plain hopping-only chain, so
# two independently-converged ground states -- one for the AUTO run, one
# for the reference run -- land on the same state rather than two
# arbitrary vectors from a degenerate manifold). "long_range=True" adds a
# single extra site-0<->site-2 hopping term on top, which is exactly the
# kind of incremental change that motivates "AUTO": the same script,
# unchanged, should keep using the cheaper "TEBD" propagator whenever the
# Hamiltonian still qualifies and fall back to "TDVP" the moment it
# doesn't, without the caller having to notice or hardcode which case
# they're in.
import time
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
from dmrgpy import timedependent

n = 6
nt, dt = 30, 0.05


def build_chain(long_range):
    fc = fermionchain.Fermionic_Chain(n)
    fc.setup_python()  # TEBD only exists on itensor_version in (3,"python")
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i]*fc.C[i+1]
    h = h + h.get_dagger()
    for i in range(n):
        h = h + 0.3*(-1)**i*fc.N[i]  # staggered field, pins a unique GS
    if long_range:
        h = h + 0.4*fc.Cdag[0]*fc.C[2] + 0.4*fc.Cdag[2]*fc.C[0]
    fc.set_hamiltonian(h)
    return fc


results = {}
for long_range in (False, True):
    fc = build_chain(long_range)
    fc.tevol_method = "AUTO"
    t0 = time.time()
    (ts, cs_auto) = timedependent.evolution_ABA(
            fc, nt=nt, dt=dt, mode="DMRG", A=fc.Cdag[0], B=fc.N[0])
    t_auto = time.time() - t0

    # cross-check against the method AUTO should have actually dispatched
    # to: "TEBD" itself for the nearest-neighbor case, "TDVP" for the one
    # that made it fall back.
    fc_ref = build_chain(long_range)
    fc_ref.tevol_method = "TDVP" if long_range else "TEBD"
    (_ts_ref, cs_ref) = timedependent.evolution_ABA(
            fc_ref, nt=nt, dt=dt, mode="DMRG", A=fc_ref.Cdag[0], B=fc_ref.N[0])

    results[long_range] = dict(ts=np.array(ts), cs_auto=np.array(cs_auto),
            cs_ref=np.array(cs_ref), t_auto=t_auto,
            ref_method=fc_ref.tevol_method)

    diff = np.max(np.abs(np.array(cs_auto) - np.array(cs_ref)))
    print(f"long_range={long_range}: AUTO matched \"{fc_ref.tevol_method}\" to "
          f"{diff:.2e} (wall {t_auto:.2f} s)")

# --- plot: AUTO's trajectory against its own reference method, for both
# the nearest-neighbor case (AUTO -> TEBD directly) and the case that
# forces the fallback (AUTO -> TDVP) ---
fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
fig.subplots_adjust(wspace=0.15, bottom=0.15)

for ax, long_range in zip(axes, (False, True)):
    r = results[long_range]
    ax.plot(r["ts"], r["cs_auto"].real, c="blue", label="AUTO")
    ax.plot(r["ts"], r["cs_ref"].real, c="red", ls="--",
            label=f'"{r["ref_method"]}" (direct)')
    title = "nearest-neighbor only (AUTO -> TEBD)" if not long_range \
            else "+ one long-range bond (AUTO -> TDVP fallback)"
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("time")
    ax.legend(fontsize=8)
axes[0].set_ylabel(r"$\langle C^\dagger_0(0)\,N_0(t)\rangle$")

plt.savefig("tevol_auto_fallback.png", dpi=150)
print("Plot saved to tevol_auto_fallback.png")
plt.show()
