# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Leaving a conserved sector *without* losing the state computed inside it
# (Many_Body_Chain.promote_to_dense; this script uses the default
# itensor_version=3 -- itensor_version="python" supports the same API,
# see examples/backend_comparison/sector_v3_VS_python).
#
# set_conserved_sector(Nf=k) confines the whole calculation to exactly k
# particles, which is what you want for the expensive part -- the sweeps.
# The price is that every operator built while it is set must itself
# conserve k: a bare c_i changes the particle number, so dmrgpy refuses to
# build it (it has to -- ITensor's AutoMPO aborts the *process* over a
# flux-violating term rather than raising anything catchable). That rules
# out precisely the quantities one usually wants from a fixed-N ground
# state: c_i|gs>, the one-body density matrix built from it, photoemission
# weights, a pairing quench.
#
# promote_to_dense() converts the block-sparse state to its exact dense
# equivalent on the chain's ordinary indices -- no re-solve, no truncation,
# nothing but the storage changes -- after which all of those are legal
# again. The sector is used for the sweeps and dropped for the
# measurement.
#
# The model is a spinless t-V chain (hopping plus nearest-neighbour
# repulsion V), whose half-filled ground state is a charge-density wave.
# Everything below is computed by actually applying c_i to a promoted
# state and checked against ED restricted to the same sector:
#   (1) |c_i|gs>|^2 = <n_i> per site -- the CDW pattern, read out through
#       an operator the sector forbids,
#   (2) the one-body density matrix <c^dag_0 c_j> vs j, from the overlaps
#       of the states c_j|gs>,
#   (3) the site-resolved quasiparticle (photoemission) weight
#       Z_i = |<gs_{N-1}| c_i |gs_N>|^2, which needs *two* sectors solved
#       separately and both promoted onto the same indices.
import numpy as np
import matplotlib.pyplot as plt

from dmrgpy import fermionchain
from dmrgpy.multioperator import MO2matrix

n = 10        # chain length
V = 2.0       # nearest-neighbour repulsion
nf = n // 2   # half filling


def tV_chain(length):
    fc = fermionchain.Fermionic_Chain(length)
    h = 0
    for i in range(length - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
        h = h + V * fc.N[i] * fc.N[i + 1]
    return fc, h, sum(fc.N)


fc, h, number = tV_chain(n)
fc.set_hamiltonian(h)
fc.maxm = 60
fc.nsweeps = 10

# --- exact reference: ED restricted to a sector ------------------------
# The number operator is diagonal in the ED product basis, so a sector is
# a set of basis states and its ground state is the lowest eigenvector of
# the corresponding submatrix, embedded back into the full space. (Do NOT
# diagonalize the full H and filter eigenvectors by <N>: eigenvalues
# degenerate across sectors come back as arbitrary superpositions.)
ed = fc.get_ED_obj()
H_ed = np.array(MO2matrix(h, ed).todense())
N_ed = np.array(MO2matrix(number, ed).todense())


def ed_sector_ground_state(target):
    sel = np.abs(np.diag(N_ed).real - target) < 1e-9
    es, vs = np.linalg.eigh(H_ed[np.ix_(sel, sel)])
    v = np.zeros(H_ed.shape[0], dtype=complex)
    v[sel] = vs[:, 0]
    return float(es[0]), v


def ed_matrix(op):
    return np.array(MO2matrix(op, ed).todense())


# --- solve at fixed N, then leave the sector ---------------------------
fc.set_conserved_sector(Nf=nf)
e_N = fc.gs_energy()
e_N_ed, v_N = ed_sector_ground_state(nf)
print("E(N=%d) = %.10f   (ED %.10f)" % (nf, e_N, e_N_ed))
assert abs(e_N - e_N_ed) < 1e-6

# c_i is illegal here, by design -- this is what the promotion is for
try:
    fc.applyoperator(fc.C[0], fc.wf0)
    raise AssertionError("a sector-violating operator was accepted")
except ValueError as err:
    print("inside the sector, c_0 is refused:\n   ", str(err).split("--")[0].strip())

fc.promote_to_dense()
print("after promote_to_dense(): conserved_sector = %s, <N> = %.10f, "
      "<H> = %.10f" % (fc.conserved_sector, fc.vev(number).real, fc.vev(h).real))
wf_N = fc.wf0.copy()

# --- (1) the CDW read out through c_i ----------------------------------
# |c_i|gs>|^2 = <gs|c^dag_i c_i|gs> = <n_i>, so applying the forbidden
# operator and taking a norm has to reproduce the site occupations.
removed = [fc.applyoperator(fc.C[i], wf_N) for i in range(n)]
weights = np.array([wf.norm() ** 2 for wf in removed])
occ_ed = np.array([float(np.conjugate(v_N) @ (ed_matrix(fc.N[i]) @ v_N)).real
                   for i in range(n)])
print("max |  |c_i|gs>|^2 - <n_i>_ED  | = %.2e" % np.max(np.abs(weights - occ_ed)))
assert np.max(np.abs(weights - occ_ed)) < 1e-6

# --- (2) the one-body density matrix from those states -----------------
# <c_i gs|c_j gs> = <gs|c^dag_i c_j|gs>, with the Jordan-Wigner strings
# carried by dmrgpy's own C/Cdag operators.
obdm = np.array([fc.overlap(removed[0], wf).real for wf in removed])
obdm_ed = np.array([float(np.conjugate(v_N) @ (ed_matrix(fc.Cdag[0] * fc.C[j])
                                               @ v_N)).real for j in range(n)])
print("max | <c^dag_0 c_j> - ED | = %.2e" % np.max(np.abs(obdm - obdm_ed)))
assert np.max(np.abs(obdm - obdm_ed)) < 1e-6

# --- (3) photoemission weight, across two sectors ----------------------
# Z_i = |<gs_{N-1}| c_i |gs_N>|^2: the probability that removing a particle
# at site i lands in the ground state of the N-1 sector rather than exciting
# it. Both states are solved in their own sector and then promoted; because
# promotion always rebases onto the *same* dense site indices (the chain's
# own, kept from construction), the two are directly comparable afterwards
# even though they were never in the same sector.
fc.set_conserved_sector(Nf=nf - 1)
e_Nm = fc.gs_energy()
e_Nm_ed, v_Nm = ed_sector_ground_state(nf - 1)
print("E(N=%d) = %.10f   (ED %.10f)" % (nf - 1, e_Nm, e_Nm_ed))
assert abs(e_Nm - e_Nm_ed) < 1e-6
fc.promote_to_dense()
wf_Nm = fc.wf0.copy()

Z = np.array([abs(fc.overlap(wf_Nm, fc.applyoperator(fc.C[i], wf_N))) ** 2
              for i in range(n)])
Z_ed = np.array([abs(np.conjugate(v_Nm) @ (ed_matrix(fc.C[i]) @ v_N)) ** 2
                 for i in range(n)])
print("max | Z_i - ED | = %.2e" % np.max(np.abs(Z - Z_ed)))
assert np.max(np.abs(Z - Z_ed)) < 1e-6
print("total removal weight sum_i Z_i = %.6f  (of sum_i <n_i> = %.6f)"
      % (Z.sum(), weights.sum()))

# --- plots -------------------------------------------------------------
sites = np.arange(n)
fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))

ax = axes[0]
ax.plot(sites, occ_ed, "o-", label=r"ED (sector-restricted) $\langle n_i\rangle$")
ax.plot(sites, weights, "x--", label=r"DMRG $\||c_i|gs\rangle\|^2$ after promotion")
ax.set_xlabel("site $i$")
ax.set_ylabel(r"$\langle n_i \rangle$")
ax.set_title("CDW read out with a sector-forbidden operator\n"
             "(n=%d, V=%g, N=%d)" % (n, V, nf))
ax.legend()

ax = axes[1]
ax.plot(sites, obdm_ed, "o-", label="ED (sector-restricted)")
ax.plot(sites, obdm, "x--", label=r"DMRG, from $\langle c_0 gs|c_j gs\rangle$")
ax.axhline(0.0, color="grey", lw=0.6)
ax.set_xlabel("site $j$")
ax.set_ylabel(r"$\langle c^\dagger_0 c_j \rangle$")
ax.set_title("One-body density matrix")
ax.legend()

ax = axes[2]
ax.plot(sites, Z_ed, "o-", label="ED (sector-restricted)")
ax.plot(sites, Z, "x--", label="DMRG, two promoted sectors")
ax.set_xlabel("site $i$")
ax.set_ylabel(r"$Z_i=|\langle gs_{N-1}|c_i|gs_N\rangle|^2$")
ax.set_title("Photoemission weight,\n$N=%d \\rightarrow N=%d$" % (nf, nf - 1))
ax.legend()

fig.tight_layout()
plt.show()
