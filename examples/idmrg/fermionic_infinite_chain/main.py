# Fermionic infinite chains (iDMRG/VUMPS) on both backends, against the
# exact free-fermion band integral.
#
# Infinite_Many_Body_Chain accepts fermionic site codes (0 = spinless,
# 1 = native spinful/Electron), and both backends now thread the
# Jordan-Wigner string themselves, locally between each term's own two
# endpoints (pyitensor/idmrg.py's _term_site_matrices,
# mpscpp3/chain_session.h's idmrg_classify_terms + idmrg_build_row).
# Before that, infinitechain.py handed the C++ backend the *finite*-chain
# Jordan-Wigner form, whose strings start at site 1 of the chain: fine as
# long as every fermionic term connected adjacent sites (the string then
# collapses to exactly the right endpoint factor with nothing in between),
# but any longer-ranged term made the backend reject the whole
# Hamiltonian.
#
# The model here is deliberately the smallest one that needs a real
# string: a dimerized (SSH-like) spinless chain with two sites A,B per
# cell,
#
#     H = t1 sum_n (A_n^dag B_n + h.c.)          intra-cell
#       + t2 sum_n (B_n^dag A_{n+1} + h.c.)      inter-cell
#       + t3 sum_n (A_n^dag A_{n+1} + h.c.)      <- skips over site B
#
# The t3 term's two endpoints have a site strictly between them, so its
# bond carries a Jordan-Wigner string (carry_ferm=True); t1/t2 alone do
# not. Being quadratic, the model's exact ground-state energy density is
# a band integral (all negative-energy single-particle levels filled),
# which is what the DMRG numbers are plotted against below.

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt

from dmrgpy import cppext
from dmrgpy import infinitechain

SPINLESS = 0   # spinless fermion site (C/Cdag/N/F)
SPINFUL = 1    # native spinful site (Cup/Cdagup/Cdn/Cdagdn/Nup/Ndn), dim 4

t1, t2 = 1.0, 0.4
MAXM = 12


def exact_density(t3, nk=20001):
    """Energy per site of the filled negative bands of the 2x2 Bloch
    matrix H(k) = [[2 t3 cos k, t1 + t2 e^{-ik}], [c.c., 0]]."""
    k = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    H = np.zeros((len(k), 2, 2), dtype=complex)
    H[:, 0, 0] = 2 * t3 * np.cos(k)
    H[:, 0, 1] = t1 + t2 * np.exp(-1j * k)
    H[:, 1, 0] = np.conj(H[:, 0, 1])
    w = np.linalg.eigvalsh(H)
    return np.where(w < 0, w, 0.0).sum(axis=1).mean() / 2


def dmrg_density(t3, itensor_version, gs_method, spinful=False):
    site = SPINFUL if spinful else SPINLESS
    ic = infinitechain.Infinite_Many_Body_Chain([site, site],
                                                itensor_version=itensor_version)
    ic.gs_method = gs_method
    ic.maxm = MAXM
    ic.maxiter = 40
    H = 0
    # A native spinful site carries both flavours at once, so the same
    # hoppings are simply written twice; the two flavours decouple, and the
    # energy density comes out exactly twice the spinless one.
    for suffix in (["up", "dn"] if spinful else [""]):
        C = [ic.get_operator("C" + suffix, i, "C") for i in range(2)]
        Cd = [ic.get_operator("Cdag" + suffix, i, "C") for i in range(2)]
        CR = [ic.get_operator("C" + suffix, i, "R") for i in range(2)]
        CdR = [ic.get_operator("Cdag" + suffix, i, "R") for i in range(2)]
        H = H + t1 * (Cd[0] * C[1] + Cd[1] * C[0])
        H = H + t2 * (Cd[1] * CR[0] + CdR[0] * C[1])
        H = H + t3 * (Cd[0] * CR[0] + CdR[0] * C[0])
    ic.set_hamiltonian(H)
    return ic.gs_energy()


t3s = np.linspace(0.0, 0.3, 7)
exact = np.array([exact_density(t3) for t3 in t3s])

runs = [("python", "idmrg"), ("python", "vumps")]
if cppext.available(3):
    runs += [(3, "idmrg"), (3, "vumps")]
else:
    print("mpscpp3 not compiled -- running the pure-Python backend only")

results = {}
for version, method in runs:
    label = "{} / {}".format(version, method)
    results[label] = np.array([dmrg_density(t3, version, method) for t3 in t3s])
    print("{:16s} max |error| vs exact = {:.2e}".format(
        label, np.abs(results[label] - exact).max()))

# The same free chain on native spinful sites: two decoupled flavours, so
# the exact answer is exactly twice the spinless one. This is what
# exercises the Electron site's own on-site spin convention.
spinful = np.array([dmrg_density(t3, "python", "idmrg", spinful=True)
                     for t3 in t3s])
print("{:16s} max |error| vs 2x exact = {:.2e}".format(
    "spinful", np.abs(spinful - 2 * exact).max()))

# The same string machinery, seen in an observable instead of an energy:
# <Cdag_0 C_r> needs a Jordan-Wigner F on every site strictly between the
# two operators. Without it the result is not "slightly off" -- it is a
# different quantity, and the discrepancy GROWS with r, which is exactly
# what corrupts a measured decay rate. The exact answer for a quadratic
# Hamiltonian is the one-body density matrix of the filled negative levels.
def exact_correlator(t3, L=200):
    N = 2 * L
    Hm = np.zeros((N, N))
    for n in range(L):
        a, b, a2 = 2 * n, 2 * n + 1, (2 * n + 2) % N
        Hm[a, b] += t1; Hm[b, a] += t1
        Hm[b, a2] += t2; Hm[a2, b] += t2
        Hm[a, a2] += t3; Hm[a2, a] += t3
    w, v = np.linalg.eigh(Hm)
    occ = v[:, w < 0]
    P = (occ.conj() @ occ.T).real
    i0 = 2 * (L // 2)
    return np.array([P[i0, i0 + r] for r in range(RMAX)])


RMAX = 9
t3_corr = 0.1
corr_exact = exact_correlator(t3_corr)
corr = {}
# itensor_version=3 with gs_method="idmrg" is the one combination with no
# correlator machinery at all (that solver's result keeps no per-sublattice
# tensor list to build one from) -- see gs_energy's own docstring. Every
# other combination supports correlators, fermionic ones included.
corr_runs = [(v, m) for v, m in runs if not (v == 3 and m == "idmrg")]
for version, method in corr_runs:
    label = "{} / {}".format(version, method)
    ic = infinitechain.Infinite_Many_Body_Chain([SPINLESS, SPINLESS],
                                                 itensor_version=version)
    ic.gs_method = method; ic.maxm = MAXM; ic.maxiter = 40
    C = [ic.get_operator("C", i, "C") for i in range(2)]
    Cd = [ic.get_operator("Cdag", i, "C") for i in range(2)]
    CR = [ic.get_operator("C", i, "R") for i in range(2)]
    CdR = [ic.get_operator("Cdag", i, "R") for i in range(2)]
    ic.set_hamiltonian(t1*(Cd[0]*C[1] + Cd[1]*C[0])
                        + t2*(Cd[1]*CR[0] + CdR[0]*C[1])
                        + t3_corr*(Cd[0]*CR[0] + CdR[0]*C[0]))
    ic.gs_energy()
    corr[label] = np.array([complex(ic.correlator("Cdag", 0, "C", r)).real
                             for r in range(RMAX)])
    print("{:16s} max |correlator error| = {:.2e}".format(
        label, np.abs(corr[label] - corr_exact).max()))

fig, (ax, axerr, axcorr) = plt.subplots(3, 1, figsize=(7, 10),
                                         gridspec_kw=dict(height_ratios=[2, 1, 2]))
ax.plot(t3s, exact, "k-", lw=2, label="exact (band integral)")
for label, e in results.items():
    ax.plot(t3s, e, "o--", ms=5, label=label)
ax.set_ylabel("energy density (per site)")
ax.set_title("Dimerized spinless chain, $t_1$={}, $t_2$={}, maxm={}".format(
    t1, t2, MAXM))
ax.legend()

for label, e in results.items():
    axerr.semilogy(t3s, np.abs(e - exact) + 1e-16, "o--", ms=5, label=label)
axerr.semilogy(t3s, np.abs(spinful - 2 * exact) + 1e-16, "s--", ms=5,
                label="spinful sites (vs 2x exact)")
axerr.set_xlabel("$t_3$ (next-cell hopping -- the term needing a JW string)")
axerr.set_ylabel("|error|")
axerr.legend(fontsize=8)

rs = np.arange(RMAX)
axcorr.plot(rs, corr_exact, "k-", lw=2, label="exact (one-body density matrix)")
for label, c in corr.items():
    axcorr.plot(rs, c, "o--", ms=5, label=label)
axcorr.axhline(0.0, color="0.8", lw=0.8, zorder=0)
axcorr.set_xlabel("separation $r$ (sites)")
axcorr.set_ylabel(r"$\langle C^\dagger_0 C_r \rangle$")
axcorr.set_title("Fermionic correlator ($t_3$={}): the Jordan-Wigner string\n"
                  "is threaded across every site between the two operators".format(t3_corr))
axcorr.legend(fontsize=8)

plt.tight_layout()
plt.show()
