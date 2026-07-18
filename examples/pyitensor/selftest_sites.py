# Self-check for dmrgpy.pyitensor Phase 2 (site type operator tables).
# There's no compiled ITensor v3 extension available to cross-check against
# directly in this environment (see cppext.available(3)), so instead this
# verifies each site's *physics*: standard operator algebra identities
# (hermiticity, commutation relations, (anti)commutators, raising/lowering
# consistency) that any correct transcription of the site headers must
# satisfy, independent of any particular reference implementation.
import numpy as np

from dmrgpy.pyitensor.sites import (SpinHalfSite, SpinOneSite, SpinThreeHalfSite,
                                     SpinTwoSite, SpinFiveHalfSite, FermionSite,
                                     ElectronSite, BosonFourSite, Z3Site, Z4Site)


def check(name, cond):
    status = "OK" if cond else "FAIL"
    print("[{}] {}".format(status, name))
    if not cond:
        raise SystemExit("selftest failed: " + name)


def std(site, opname):
    """The standard (physicist) matrix for an operator: the transpose of
    this engine's own (in,out) storage convention -- see sites/base.py's
    module docstring for why the transpose is the right conversion."""
    return site.matrix(opname).T


def is_hermitian(M, tol=1e-10):
    return np.allclose(M, M.conj().T, atol=tol)


def commutator(A, B):
    return A @ B - B @ A


def anticommutator(A, B):
    return A @ B + B @ A


def test_spin_site(site, spin, has_pm=True, has_s2_op=True, name=""):
    Sx, Sy, Sz = std(site, "Sx"), std(site, "Sy"), std(site, "Sz")
    check(name + ": Sx hermitian", is_hermitian(Sx))
    check(name + ": Sy hermitian", is_hermitian(Sy))
    check(name + ": Sz hermitian", is_hermitian(Sz))
    check(name + ": [Sx,Sy]=iSz", np.allclose(commutator(Sx, Sy), 1j * Sz))
    check(name + ": [Sy,Sz]=iSx", np.allclose(commutator(Sy, Sz), 1j * Sx))
    check(name + ": [Sz,Sx]=iSy", np.allclose(commutator(Sz, Sx), 1j * Sy))
    S2 = Sx @ Sx + Sy @ Sy + Sz @ Sz
    check(name + ": Sx^2+Sy^2+Sz^2 = s(s+1)*I", np.allclose(
        S2, spin * (spin + 1) * np.eye(site.dim)))
    if has_pm:
        Sp, Sm = std(site, "Sp"), std(site, "Sm")
        check(name + ": Sp = Sx + i*Sy", np.allclose(Sp, Sx + 1j * Sy))
        check(name + ": Sm = Sx - i*Sy", np.allclose(Sm, Sx - 1j * Sy))
        check(name + ": S+ alias matches Sp", np.allclose(std(site, "S+"), Sp))
        check(name + ": S- alias matches Sm", np.allclose(std(site, "S-"), Sm))


def test_fermion_site():
    C, Cdag, N, F = std(FermionSite, "C"), std(FermionSite, "Cdag"), std(FermionSite, "N"), std(FermionSite, "F")
    check("Fermion: Cdag = C^dagger", np.allclose(Cdag, C.conj().T))
    check("Fermion: N = Cdag C", np.allclose(N, Cdag @ C))
    check("Fermion: {C,Cdag} = Id", np.allclose(anticommutator(C, Cdag), np.eye(2)))
    check("Fermion: C^2 = 0", np.allclose(C @ C, 0))
    check("Fermion: F = Id - 2N", np.allclose(F, np.eye(2) - 2 * N))
    check("Fermion: A/Adag alias C/Cdag", np.allclose(std(FermionSite, "A"), C)
          and np.allclose(std(FermionSite, "Adag"), Cdag))


def test_electron_site():
    s = ElectronSite
    Cup, Cdagup = std(s, "Cup"), std(s, "Cdagup")
    Cdn, Cdagdn = std(s, "Cdn"), std(s, "Cdagdn")
    Nup, Ndn, Ntot = std(s, "Nup"), std(s, "Ndn"), std(s, "Ntot")
    F, Fup, Fdn = std(s, "F"), std(s, "Fup"), std(s, "Fdn")
    Id4 = np.eye(4)

    check("Electron: Cdagup = Cup^dagger", np.allclose(Cdagup, Cup.conj().T))
    check("Electron: Cdagdn = Cdn^dagger", np.allclose(Cdagdn, Cdn.conj().T))
    check("Electron: Nup = Cdagup Cup", np.allclose(Nup, Cdagup @ Cup))
    check("Electron: Ndn = Cdagdn Cdn", np.allclose(Ndn, Cdagdn @ Cdn))
    check("Electron: Ntot = Nup + Ndn", np.allclose(Ntot, Nup + Ndn))
    check("Electron: {Cup,Cdagup} = Id", np.allclose(anticommutator(Cup, Cdagup), Id4))
    check("Electron: {Cdn,Cdagdn} = Id", np.allclose(anticommutator(Cdn, Cdagdn), Id4))
    check("Electron: {Cup,Cdn} = 0 (same-site fermion ordering sign)",
          np.allclose(anticommutator(Cup, Cdn), 0))
    check("Electron: Cup^2 = 0", np.allclose(Cup @ Cup, 0))
    check("Electron: Cdn^2 = 0", np.allclose(Cdn @ Cdn, 0))
    # F = (-1)^Ntot, not Id-2*Ntot -- that only agrees for Ntot in {0,1};
    # the doubly-occupied state (Ntot=2) needs the periodic form.
    check("Electron: F = (-1)^Ntot", np.allclose(F, np.diag((-1.0) ** np.diag(Ntot))))
    check("Electron: F = Fup*Fdn (both diagonal)", np.allclose(F, Fup @ Fdn))

    Sz, Sp, Sm = std(s, "Sz"), std(s, "Sp"), std(s, "Sm")
    check("Electron: Sz = (Nup-Ndn)/2", np.allclose(Sz, 0.5 * (Nup - Ndn)))
    check("Electron: Sp = Cdagup Cdn", np.allclose(Sp, Cdagup @ Cdn))
    check("Electron: Sm = Cdagdn Cup", np.allclose(Sm, Cdagdn @ Cup))


def test_boson_site():
    A, Adag, N = std(BosonFourSite, "A"), std(BosonFourSite, "Adag"), std(BosonFourSite, "N")
    check("Boson4: Adag = A^dagger", np.allclose(Adag, A.conj().T))
    check("Boson4: N = Adag A", np.allclose(N, Adag @ A))
    check("Boson4: N eigenvalues are 0,1,2,3", np.allclose(np.diag(N), [0, 1, 2, 3]))
    # [a,adag] = 1 everywhere except the truncated top level (level 3 has no
    # level 4 to raise into, so the usual boson commutator only holds on the
    # first dim-1 levels of a level-truncated Fock space).
    comm = commutator(A, Adag)
    check("Boson4: [A,Adag]=Id on the untruncated levels", np.allclose(np.diag(comm)[:3], 1.0))


def test_zn_site(site, n, name):
    Sig, SigDag = std(site, "Sig"), std(site, "SigDag")
    Tau, TauDag = std(site, "Tau"), std(site, "TauDag")
    Id = np.eye(site.dim)
    check(name + ": SigDag = Sig^dagger", np.allclose(SigDag, Sig.conj().T))
    check(name + ": Sig is unitary (Sig SigDag = Id)", np.allclose(Sig @ SigDag, Id))
    power = np.eye(site.dim, dtype=complex)
    for _ in range(n):
        power = power @ Sig
    check(name + ": Sig^{} = Id (order-{} shift/clock op)".format(n, n), np.allclose(power, Id))
    check(name + ": TauDag = Tau^dagger", np.allclose(TauDag, Tau.conj().T))
    check(name + ": Tau is unitary", np.allclose(Tau @ TauDag, Id))
    omega = np.exp(2j * np.pi / n)
    check(name + ": Sig*Tau = omega*Tau*Sig (clock algebra)",
          np.allclose(Sig @ Tau, omega * (Tau @ Sig)))


if __name__ == "__main__":
    test_spin_site(SpinHalfSite, 0.5, name="SpinHalf")
    test_spin_site(SpinOneSite, 1.0, name="SpinOne")
    test_spin_site(SpinThreeHalfSite, 1.5, has_pm=False, name="SpinThreeHalf")
    test_spin_site(SpinFiveHalfSite, 2.5, has_pm=False, name="SpinFiveHalf")
    test_spin_site(SpinTwoSite, 2.0, name="SpinTwo")
    test_fermion_site()
    test_electron_site()
    test_boson_site()
    test_zn_site(Z3Site, 3, "Z3")
    test_zn_site(Z4Site, 4, "Z4")
    print("All pyitensor Phase 2 checks passed.")
