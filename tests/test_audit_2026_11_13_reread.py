"""Regressions of the 2026-11-13 Opus re-read of the Kondo commits of
2026-11-10 to 13: the SU(2) reduction of submode="KRYLOV" on a chain of
spinful fermions with both kinds of total spin, the n_gs manifold of the
TASK route at a temperature below the truncation split of its members, and
the memory budget of order="LL" through the public call."""
import warnings

import numpy as np
import pytest

from dmrgpy import fermionchain, spinchain
from dmrgpy.kondospectrumtk import conductance, edkondo, krylovkondo, leadinglog

from test_kondo_task_convergence import VERSIONS, _odd_chain


def _hubbard(n=3):
    """A Hubbard chain with no Nf sector and no field: even and odd
    electron numbers, so integer and half-integer total spins together."""
    fc = fermionchain.Spinful_Fermionic_Chain_Native(n)
    h = 0
    for i in range(n - 1):
        h = h - (fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1])
    h = h + h.get_dagger()
    for i in range(n):
        h = h + 8.*(fc.Nup[i] - 0.5)*(fc.Ndn[i] - 0.5)
    fc.set_hamiltonian(h)
    return fc


def test_su2_keeps_the_multiplets_of_both_kinds_of_spin():
    """The smallest |Sz_tot| is 0 here, and taking only those states
    dropped every half-integer multiplet: 0.338 of the peak off ED."""
    fc = _hubbard()
    eV = np.linspace(-0.6, 0.6, 9) + 1e-3
    kw = dict(site=0, Jrho_s=-0.05, U=0.2, T=0.3, kB=1., omega0=4., Gamma0=1e-2)
    _, ed = fc.get_kondo_spectrum(eV, **kw)
    _, kr = fc.get_kondo_spectrum(eV, submode="KRYLOV", nkry=40, p_cut=1e-14,
                                  su2=True, **kw)
    assert np.max(np.abs(kr - ed)) < 1e-11*np.max(np.abs(ed))
    spins = krylovkondo.KrylovKondoSpectrum(fc, 0, 0.3, kB=1., su2=True).spins
    assert np.any(np.mod(2*spins, 2) == 0) and np.any(np.mod(2*spins, 2) == 1)


@pytest.mark.parametrize("version", VERSIONS)
def test_a_truncated_doublet_is_averaged_below_its_split(version):
    """At 1 mK, 1e-3 kT is below the 1e-9 eV at which truncation puts the
    partner of a member, and the call was refused while T=0 and 0.1 K ran;
    the partner is a member within Gamma0, the split the n_gs average
    accepts."""
    sc = _odd_chain(version)
    eV = np.linspace(-4e-3, 4e-3, 41)
    kw = dict(site=0, Jrho_s=-0.1, U=0.2, T=1e-3, omega0=20e-3, Gamma0=5e-5,
              order=3)
    _, ed = sc.get_kondo_spectrum(eV, mode="ED", **kw)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # the split of the truncated doublet
        _, tk = sc.get_kondo_spectrum(eV, mode="DMRG", nkry=60, n_gs=2, **kw)
        # one member alone is still not the thermal state
        with pytest.raises(ValueError):
            sc.get_kondo_spectrum(eV, mode="DMRG", nkry=60, **kw)
    assert np.max(np.abs(tk - ed)) < 1e-3*np.max(np.abs(ed))


def test_mem_gb_reaches_the_leading_log_flow(monkeypatch):
    """order="LL" at T>0 held every bias in one flow on the CPU, and the
    public call had no way to pass mem_gb."""
    sc = spinchain.Spin_Chain(["S=1/2"]*3)
    h = 0
    for i in range(2):
        h = h + 1e-3*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    for i in range(3):
        h = h + 0.3e-3*sc.Sz[i]
    sc.set_hamiltonian(h)
    eV = np.linspace(-3e-3, 3e-3, 24) + 1e-5
    kw = dict(site=0, Jrho_s=-0.05, U=0.2, T=2.0, order="LL")
    ks = edkondo.KondoSpectrum(sc, 0, 2.0)
    fixed, per_w = leadinglog.memory_model(len(ks.e), len(conductance._occupied_states(ks)))
    mem_gb = (fixed + 16*per_w)*1.01/1e9  # 8 biases per flow, so 3 flows
    chunks = []
    original = leadinglog.leading_log_dIdV
    def spy(*a, **k):
        out, info = original(*a, return_info=True, **k)
        chunks.append(info["n_chunks"])
        return out
    monkeypatch.setattr(leadinglog, "leading_log_dIdV", spy)
    _, one = sc.get_kondo_spectrum(eV, **kw)
    _, cut = sc.get_kondo_spectrum(eV, mem_gb=mem_gb, **kw)
    assert chunks == [1, 3]
    assert np.max(np.abs(cut/one - 1)) < 1e-9
    with pytest.raises(ValueError, match="mem_gb"):
        sc.get_kondo_spectrum(eV, mem_gb=1e-9, **kw)
    with pytest.raises(TypeError, match="mem_gb"):
        sc.get_kondo_spectrum(eV, site=0, Jrho_s=-0.05, T=2.0, order=3, mem_gb=1.)
