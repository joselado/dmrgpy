"""The native correlation matrix <Cdag_i C_j> on itensor_version=
"julia_live" (dmmode="full", mpsjulialive/correlations.jl) against
dmmode="explicit" on the same state, "python" dmmode="full" and ED. The
chain has a next-nearest-neighbor hopping and a staggered onsite energy,
so the odd-distance entries are not zero by symmetry and a wrong
Jordan-Wigner sign would show in them."""
import types

import numpy as np
import pytest

from dmrgpy import fermionchain, spinchain
from dmrgpy.entropytk import correlationentropy

from _helpers import julia_live_param, setup_backend


def chain(version, n=8):
    fc = fermionchain.Fermionic_Chain(n)
    if version == "ED": fc.mode = "ED"
    else: setup_backend(fc, version)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
        h = h + 0.5 * (fc.N[i] - 0.5) * (fc.N[i + 1] - 0.5)
    for i in range(n - 2):
        h = h + 0.4 * (fc.Cdag[i] * fc.C[i + 2] + fc.Cdag[i + 2] * fc.C[i])
    for i in range(n):
        h = h + (0.3 * (-1)**i + 0.1) * fc.N[i]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps = 40, 20
    return fc


@pytest.fixture(scope="module")
def reference():
    ed = chain("ED")
    return np.asarray(ed.get_correlation_matrix(dmmode="explicit"))


def test_reference_has_odd_distance_entries(reference):
    odd = [reference[i, i + 1] for i in range(7)] + [reference[0, 3]]
    assert np.min(np.abs(odd)) > 1e-3


@pytest.mark.parametrize("version", ["python", julia_live_param()])
def test_full_matches_ed_and_explicit(version, reference):
    fc = chain(version)
    wf = fc.get_gs()
    full = np.asarray(fc.get_correlation_matrix(dmmode="full", wf=wf))
    explicit = np.asarray(fc.get_correlation_matrix(dmmode="explicit", wf=wf))
    assert np.max(np.abs(full - explicit)) < 1e-8
    assert np.max(np.abs(full - reference)) < 1e-8
    # an unnormalized input gives the matrix of the normalized state
    full2 = np.asarray(fc.get_correlation_matrix(dmmode="full", wf=2.0 * wf))
    assert np.max(np.abs(full2 - full)) < 1e-10


@pytest.mark.parametrize("version", [julia_live_param()])
def test_julia_default_under_a_sector_is_full(version):
    """_default_dmmode picks "full" for a julia_live state of spinless
    fermions when a sector is set, and not for Nambu or a spin chain."""
    fc = chain(version)
    wf = fc.get_gs()
    sector = types.SimpleNamespace(conserved_sector={"Nf": 4})
    assert correlationentropy._default_dmmode(sector, wf, "electron") == "full"
    assert correlationentropy._default_dmmode(sector, wf, "Nambu") == "explicit"
    sc = spinchain.Spin_Chain(["S=1/2"] * 4)
    setup_backend(sc, version)
    sc.set_hamiltonian(sc.Sz[0] * sc.Sz[1] + sc.Sx[1] * sc.Sx[2])
    assert correlationentropy._default_dmmode(sector, sc.get_gs(),
                                              "electron") == "explicit"
