"""Regressions for docs/audit_2026_11_21_hole_hunt.md, findings 6 and 7
(cluster metts): the METTS collapse basis guard runs on every backend and
in both METTS calls, and it refuses every basis whose samples stay in one
sector of a charge or a parity that H conserves."""
import numpy as np
import pytest

from dmrgpy import cppext, fermionchain, spinchain

from _helpers import julia_live_param

BACKENDS = [
    "python",
    pytest.param(3, marks=pytest.mark.skipif(
        not cppext.available(3),
        reason="requires the compiled mpscpp3 (ITensor v3) extension")),
    julia_live_param(),
]

KW = dict(maxm=16, nsweeps=6, cutoff=1e-12)
T = 0.5


@pytest.fixture(autouse=True)
def _seed():
    np.random.seed(0) # the random witness of the guard


def heisenberg(v, n=4):
    """conserves Sxtot, Sytot and Sztot"""
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=v, **KW)
    h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
            for i in range(n-1))
    sc.set_hamiltonian(h)
    return sc, h


def ising_x(v, n=4):
    """transverse-field Ising chain with its field along x: conserves
    Px = prod_i 2Sx_i and no sum of one-site operators"""
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=v, **KW)
    h = sum(sc.Sz[i]*sc.Sz[i+1] for i in range(n-1)) + 0.5*sum(sc.Sx)
    sc.set_hamiltonian(h)
    return sc, h


def hopping(v, n=4):
    """conserves N"""
    fc = fermionchain.Fermionic_Chain(n, itensor_version=v, **KW)
    h = sum(-1.0*(fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i])
            for i in range(n-1))
    h = h + 0.3*fc.N[0] - 0.2*fc.N[n-1] + 0.1*fc.N[1]
    fc.set_hamiltonian(h)
    return fc, h


def kitaev(v, n=4):
    """a pairing term: conserves (-1)^N and not N"""
    fc = fermionchain.Fermionic_Chain(n, itensor_version=v, **KW)
    h = sum(-1.0*(fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i])
            + 0.6*(fc.C[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.Cdag[i])
            for i in range(n-1))
    h = h + 0.3*sum(fc.N)
    fc.set_hamiltonian(h)
    return fc, h


# (chain, basis, what the message names). Each basis gives product states
# of definite charge or parity, so before the fix the samples stayed in the
# sector of the first one, with stderr 0 on that charge.
FROZEN = [
    (heisenberg, ("Sz",), "conserves"),
    (heisenberg, ("Sx",), "conserves"), # not diagonal, Sxtot conserved
    (ising_x, ("Sx",), "the parity prod_i 2Sx_i"),
    (hopping, ("N",), "conserves"),
    (hopping, ("n",), "conserves"), # the diagonal aliases of N
    (hopping, ("projOcc",), "conserves"),
    (hopping, ("projEmp",), "conserves"),
    (hopping, ("FermiPhase",), "conserves"),
    (kitaev, ("N",), "fermion parity"),
    (heisenberg, ("S-",), "not Hermitian"),
]
FROZEN_IDS = ["heis-Sz", "heis-Sx", "isingx-Sx", "hop-N", "hop-n",
              "hop-projOcc", "hop-projEmp", "hop-FermiPhase", "kitaev-N",
              "heis-S-"]


# Finding 6 (every backend, both calls) and finding 7 (one non-diagonal
# name, the Z2 parities, the aliases of N)

@pytest.mark.parametrize("chain,basis,match", FROZEN, ids=FROZEN_IDS)
@pytest.mark.parametrize("v", BACKENDS)
def test_metts_vev_refuses_a_frozen_basis(v, chain, basis, match):
    sc, h = chain(v)
    with pytest.raises(ValueError, match=match) as err:
        sc.metts_vev(h, T, nsamples=4, nwarmup=1, basis_ops=basis, seed=1)
    assert str(err.value).startswith("metts_vev:")


@pytest.mark.parametrize("chain,basis,match", FROZEN, ids=FROZEN_IDS)
@pytest.mark.parametrize("v", BACKENDS)
def test_metts_dynamical_correlator_refuses_a_frozen_basis(v, chain, basis,
                                                           match):
    sc, h = chain(v)
    with pytest.raises(ValueError, match=match) as err:
        sc.metts_dynamical_correlator((h, h), T, nt=2, dt=0.1, nsamples=2,
                                      nwarmup=1, basis_ops=basis, seed=1)
    assert str(err.value).startswith("metts_dynamical_correlator:")


@pytest.mark.parametrize("v", BACKENDS)
def test_a_mixing_basis_passes_the_guard(v):
    """The default basis, and ("Sz",) on the Ising chain with its field
    along x, where neither Sztot nor prod_i 2Sz_i is conserved (the
    control of the 2026-10-23 guard)"""
    from dmrgpy.vevtk.mettsvev import _check_basis
    sc, h = heisenberg(v)
    _check_basis(sc, ("Sz", "Sx"))
    _check_basis(sc, ("Sz", "Sx"), caller="metts_dynamical_correlator")
    sc, h = ising_x(v)
    _check_basis(sc, ("Sz", "Sx"))
    _check_basis(sc, ("Sz",))


@pytest.mark.parametrize("v", BACKENDS)
def test_a_mixing_basis_samples_every_sector_and_matches_ed(v):
    """The default basis ("Sz","Sx") on the Heisenberg chain of the
    record's controls (4 sites, T=0.5, 10 samples, seed 1): <H> within 4
    error bars of ED (the record has 2.1, 0.7 and 0.4 of them on "python",
    v3 and julia_live, for the same chain and seed), and a finite error
    bar on both Sztot^2 and Sxtot^2, where a frozen chain gives 1e-15 on
    its charge"""
    sc, h = heisenberg(v)
    sz, sx = sum(sc.Sz), sum(sc.Sx)
    (_, ez), (_, ex), (e, ee) = sc.metts_vev([sz*sz, sx*sx, h], T,
            nsamples=10, nwarmup=2, seed=1)
    e_ed = sc.vev(h, mode="ED", T=T).real
    assert abs(e.real - e_ed) < 4*ee
    assert ez > 1e-3 and ex > 1e-3


# Finding 6, defence in depth: the eigenbasis of the collapse operator

def test_python_eigenbasis_refuses_a_non_hermitian_name():
    from dmrgpy.pyitensor import metts
    from dmrgpy.pyitensor.sites import SiteX
    sites = SiteX([2, 2])
    evals, _ = metts._eigenbasis(sites, "Sx", 1)
    assert np.allclose(evals, [-0.5, 0.5])
    with pytest.raises(ValueError, match="not Hermitian"):
        metts._eigenbasis(sites, "S-", 1)


@pytest.mark.parametrize("v", [julia_live_param()])
def test_julia_eigenbasis_refuses_a_non_hermitian_name(v):
    """eigen(Hermitian(M)) read the upper triangle of S-, which is zero"""
    from dmrgpy.mpsjulialive.juliasession import Main, to_julia_strvec
    sc = spinchain.Spin_Chain(["S=1/2"]*3, itensor_version=v)
    Main.metts_build_eigcache(sc.jlsites, to_julia_strvec(["Sx"]))
    with pytest.raises(Exception, match="not Hermitian"):
        Main.metts_build_eigcache(sc.jlsites, to_julia_strvec(["S-"]))
