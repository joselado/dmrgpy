# Self-check for dmrgpy.pyitensor Phase 3 (AutoMPO term compiler / Jordan-
# Wigner threading). Direct matrix comparison against dmrgpy's independent
# ED backend (pyfermion.MBFermion) isn't meaningful here since the two use
# different basis-ordering conventions for the same occupation-number
# Hilbert space (MBFermion's states.generate_basis is little-endian in
# site index, this module's dense_matrix() is big-endian, see the
# discussion this test's docstrings would otherwise duplicate) -- so
# instead this compares *eigenvalues* of matching Hamiltonians, which is
# basis-ordering-independent and just as sensitive to a Jordan-Wigner sign
# bug (a wrong sign generically changes the spectrum).
import numpy as np

from dmrgpy.pyitensor.autompo import AutoMPO
from dmrgpy.pyitensor.sites import SiteX
from dmrgpy.pyfermion.mbfermion import MBFermion


def check(name, cond):
    status = "OK" if cond else "FAIL"
    print("[{}] {}".format(status, name))
    if not cond:
        raise SystemExit("selftest failed: " + name)


def eigs(M):
    return np.sort(np.linalg.eigvalsh(M))


def test_spinless_fermion_jw():
    n = 5
    sites = SiteX([0] * n)  # 0 = spinless fermion
    ampo = AutoMPO(sites)
    # nearest-neighbor hopping
    for i in range(1, n):
        ampo.add(1.0, "Cdag", i, "C", i + 1)
        ampo.add(1.0, "Cdag", i + 1, "C", i)
    # a range-2 hop (skips a site -- exercises a 1-site-long F string)
    ampo.add(0.5, "Cdag", 2, "C", 4)
    ampo.add(0.5, "Cdag", 4, "C", 2)
    # a long-range hop spanning the whole chain (exercises a long F string)
    ampo.add(0.3, "Cdag", 1, "C", 5)
    ampo.add(0.3, "Cdag", 5, "C", 1)
    # an explicit 4-fermion-operator term, mirroring
    # chain_session.h::four_correlation_tensor's Cdag_i C_j Cdag_k C_l shape
    # (with its own +h.c.) -- the general N-fermion JW case, not just pairs.
    ampo.add(0.2, "Cdag", 1, "C", 2, "Cdag", 3, "C", 4)
    ampo.add(0.2, "Cdag", 4, "C", 3, "Cdag", 2, "C", 1)
    mine = ampo.dense_matrix()
    check("spinless-fermion Hamiltonian is Hermitian", np.allclose(mine, mine.conj().T))

    mb = MBFermion(n)
    ref = mb.get_zero()
    for i in range(1, n):
        ref = ref + mb.Cdag[i - 1] @ mb.C[i] + mb.Cdag[i] @ mb.C[i - 1]
    ref = ref + 0.5 * (mb.Cdag[1] @ mb.C[3] + mb.Cdag[3] @ mb.C[1])
    ref = ref + 0.3 * (mb.Cdag[0] @ mb.C[4] + mb.Cdag[4] @ mb.C[0])
    ref = ref + 0.2 * (mb.Cdag[0] @ mb.C[1] @ mb.Cdag[2] @ mb.C[3]
                        + mb.Cdag[3] @ mb.C[2] @ mb.Cdag[1] @ mb.C[0])
    ref = ref.toarray()
    check("reference (independent ED backend) Hamiltonian is Hermitian", np.allclose(ref, ref.conj().T))

    check("nearest+range-2+long-range+4-fermion-term spectrum matches "
          "dmrgpy's independent ED backend (Jordan-Wigner sign is correct)",
          np.allclose(eigs(mine), eigs(ref), atol=1e-9))


def test_spinless_fermion_single_site_number():
    # A term with no fermionic operators at all that still needs "N" to be
    # built correctly through AutoMPO's coefficient bookkeeping.
    n = 3
    sites = SiteX([0] * n)
    ampo = AutoMPO(sites)
    for i in range(1, n + 1):
        ampo.add(1.0, "N", i)
    mine = eigs(ampo.dense_matrix())

    mb = MBFermion(n)
    ref = mb.get_zero()
    for i in range(n):
        ref = ref + mb.N[i]
    check("sum of N_i has integer eigenvalues 0..n", np.allclose(mine, eigs(ref.toarray())))
    check("sum of N_i eigenvalues are exactly 0,1,1,1,2,2,2,3",
          np.allclose(sorted(mine), [0, 1, 1, 1, 2, 2, 2, 3]))


def test_hubbard_jw():
    # ElectronSite packs spin-up and spin-down into one 4-dim site, with an
    # intrinsic same-site sign baked into Cdn/Cdagdn (see sites/fermion.py's
    # docstring). Cross-check that sign, plus ordinary cross-site JW
    # threading for two independent fermion species at once, against an
    # independent flat-orbital ED reference: treat site i's up/down as two
    # separate spinless-fermion orbitals 2i and 2i+1 (0-based) in one big
    # MBFermion(2*n) -- if ElectronSite's built-in sign or this module's
    # cross-site JW threading were wrong, this Hubbard spectrum would not
    # match a model built with completely independent fermionic orbitals.
    n = 3
    sites = SiteX([1] * n)  # 1 = spinful fermion (Hubbard)
    ampo = AutoMPO(sites)
    t, U, J = 1.0, 2.0, 0.3
    for i in range(1, n):
        ampo.add(t, "Cdagup", i, "Cup", i + 1)
        ampo.add(t, "Cdagup", i + 1, "Cup", i)
        ampo.add(t, "Cdagdn", i, "Cdn", i + 1)
        ampo.add(t, "Cdagdn", i + 1, "Cdn", i)
    for i in range(1, n + 1):
        ampo.add(U, "Nup", i, "Ndn", i)
    # a spin-flip-like term mixing species across sites, for good measure
    ampo.add(J, "Cdagup", 1, "Cdn", 3)
    ampo.add(J, "Cdagdn", 3, "Cup", 1)
    mine = eigs(ampo.dense_matrix())

    mb = MBFermion(2 * n)

    def up(i):
        return 2 * (i - 1)  # 0-based flat-orbital index

    def dn(i):
        return 2 * (i - 1) + 1

    ref = mb.get_zero()
    for i in range(1, n):
        ref = ref + t * (mb.Cdag[up(i)] @ mb.C[up(i + 1)] + mb.Cdag[up(i + 1)] @ mb.C[up(i)])
        ref = ref + t * (mb.Cdag[dn(i)] @ mb.C[dn(i + 1)] + mb.Cdag[dn(i + 1)] @ mb.C[dn(i)])
    for i in range(1, n + 1):
        ref = ref + U * (mb.N[up(i)] @ mb.N[dn(i)])
    ref = ref + J * (mb.Cdag[up(1)] @ mb.C[dn(3)] + mb.Cdag[dn(3)] @ mb.C[up(1)])
    refe = eigs(ref.toarray())

    check("Hubbard chain (hopping + U + spin-flip) spectrum matches an "
          "independent flat-orbital ED reference (ElectronSite's built-in "
          "same-site sign and cross-site JW threading both correct)",
          np.allclose(mine, refe, atol=1e-9))


def test_spin_chain_no_jw_needed():
    # Spin operators are never fermionic (is_fermionic("Sz") etc. is always
    # False), so this exercises AutoMPO's term/coefficient bookkeeping with
    # Jordan-Wigner threading never triggering -- cross-checked here against
    # an explicitly hand-built Kronecker-product reference using the exact
    # same site-ordering convention (no basis-reordering ambiguity, so this
    # one *is* a direct matrix comparison).
    from dmrgpy.pyitensor.sites import SpinHalfSite
    n = 4
    sites = SiteX([2] * n)  # 2 = spin-1/2
    ampo = AutoMPO(sites)
    for i in range(1, n):
        ampo.add(1.0, "Sz", i, "Sz", i + 1)
        ampo.add(0.5, "Sp", i, "Sm", i + 1)
        ampo.add(0.5, "Sm", i, "Sp", i + 1)
    mine = ampo.dense_matrix()

    Sz, Sp, Sm, Id = (SpinHalfSite.matrix(nm).T for nm in ("Sz", "Sp", "Sm", "Id"))

    def kron_at(ops_by_site):
        full = ops_by_site.get(1, Id)
        for i in range(2, n + 1):
            full = np.kron(full, ops_by_site.get(i, Id))
        return full

    ref = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for i in range(1, n):
        ref += kron_at({i: Sz, i + 1: Sz})
        ref += 0.5 * kron_at({i: Sp, i + 1: Sm})
        ref += 0.5 * kron_at({i: Sm, i + 1: Sp})

    check("Heisenberg chain matches a hand-built Kronecker reference exactly",
          np.allclose(mine, ref))
    check("Heisenberg chain is Hermitian", np.allclose(mine, mine.conj().T))


if __name__ == "__main__":
    test_spinless_fermion_jw()
    test_spinless_fermion_single_site_number()
    test_hubbard_jw()
    test_spin_chain_no_jw_needed()
    print("All pyitensor Phase 3 checks passed.")
