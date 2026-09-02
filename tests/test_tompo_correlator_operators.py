"""`name=` accepts the operators `toMPO()` returns, not just MultiOperators.

`Many_Body_Chain.toMPO()` is the documented build-once/reuse-many fast
path for a dynamical correlator: build the two operators once, hand them
to `get_dynamical_correlator(name=(A,B))`, sweep the frequency grid.
Centralizing `name=` resolution in
`Many_Body_Chain.get_dynamical_correlator` (commit 1b87543) routed every
submode through `operatornames.str2MO`, whose explicit-pair branch only
ever admitted `MultiOperator` -- so the paths that used to consume an
already-built operator directly (`mode="ED"` and DMRG `submode="EX"`,
neither of which reached `str2MO` before) started rejecting it with a
`ValueError`, including the product form real code builds
(`toMPO(Cdag[0])*toMPO(C[1])`).

What is pinned here: the compiled operators pass through and give the
same numbers as the MultiOperators they were built from, and the
submodes that genuinely cannot take one -- they rebuild the operator
inside the backend from `MultiOperator.to_terms()` -- say so with a
`TypeError` naming `toMPO`, instead of an `AttributeError` several
frames deep.
"""
import numpy as np
import pytest

from dmrgpy import fermionchain, spinchain

ES = np.linspace(-1.0, 4.0, 40)
DELTA = 0.1


def _heisenberg(n=6):
    sc = spinchain.Spin_Chain(["S=1/2"] * n)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
                + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h)
    return sc


def _correlator(sc, A, B, es=ES, **kwargs):
    return sc.get_dynamical_correlator(name=(A, B), es=es, delta=DELTA,
                                       **kwargs)[1]


@pytest.mark.parametrize("submode", ["ED", "KPM", "INV", "CVM", "EX", "TD"])
def test_ed_accepts_tompo_operators(submode):
    """Every ED submode: toMPO(mode="ED") output is accepted and gives the
    same spectrum as the MultiOperator it was built from. ED is exact and
    deterministic here, so the two agree to machine precision rather than
    to a DMRG tolerance."""
    sc = _heisenberg()
    ref = _correlator(sc, sc.Sz[0], sc.Sz[0], mode="ED", submode=submode)
    sc = _heisenberg()
    B = sc.toMPO(sc.Sz[0], mode="ED")
    out = _correlator(sc, B, B, mode="ED", submode=submode)
    assert np.max(np.abs(np.array(out) - np.array(ref))) == pytest.approx(
        0.0, abs=1e-10)


@pytest.mark.parametrize("submode", ["CVM", "EX", "ROOTN"])
def test_dmrg_accepts_tompo_operators(submode):
    """The DMRG submodes that consume an operator by *applying* it accept
    a StaticOperator. itensor_version="python" so this needs no compiled
    extension."""
    # ROOTN costs N=8 sequential fractional-power Lanczos solves per
    # frequency, so it gets a coarse grid; the comparison is unaffected
    es = np.linspace(-1.0, 4.0, 4) if submode == "ROOTN" else ES
    sc = _heisenberg()
    sc.setup_python()
    ref = _correlator(sc, sc.Sz[0], sc.Sz[0], es=es, submode=submode)
    sc = _heisenberg()
    sc.setup_python()
    B = sc.toMPO(sc.Sz[0])
    out = _correlator(sc, B, B, es=es, submode=submode)
    # EX rediagonalizes from a fresh random start each time, so this is a
    # "same spectrum", not a "same bits", comparison
    assert np.max(np.abs(np.array(out) - np.array(ref))) < 0.1


@pytest.mark.parametrize("submode", ["KPM", "TD"])
def test_symbolic_only_submodes_name_the_problem(submode):
    """KPM/TD rebuild their operators inside the backend from to_terms(),
    so they cannot take a StaticOperator -- and must say that."""
    sc = _heisenberg()
    sc.setup_python()
    B = sc.toMPO(sc.Sz[0])
    with pytest.raises(TypeError) as excinfo:
        _correlator(sc, B, B, submode=submode)
    message = str(excinfo.value)
    assert "toMPO" in message
    assert submode in message


def test_product_of_tompo_operators():
    """The form real code builds: a product of two already-built
    operators, with its dagger as the second half of the pair (the
    cotunneling channel of a spinful-fermion forward model)."""
    def build():
        fc = fermionchain.Spinful_Fermionic_Chain(3)
        h = 0
        for i in range(2):
            h = h + fc.Cdagup[i] * fc.Cup[i + 1] + fc.Cdagdn[i] * fc.Cdn[i + 1]
        h = h + h.get_dagger()
        for i in range(3):
            h = h + 2.0 * fc.Nup[i] * fc.Ndn[i]
        fc.set_hamiltonian(h)
        return fc

    fc = build()
    Amo = fc.Cdagup[0] * fc.Cup[1]
    ref = _correlator(fc, Amo, Amo.get_dagger(), mode="ED", submode="ED")

    fc = build()
    Cdagup = [fc.toMPO(A, mode="ED") for A in fc.Cdagup]
    Cup = [fc.toMPO(A, mode="ED") for A in fc.Cup]
    A01 = Cdagup[0] * Cup[1]
    out = _correlator(fc, A01, A01.get_dagger(), mode="ED", submode="ED")

    assert np.max(np.abs(np.array(ref))) > 1e-3  # the channel is not empty
    assert np.max(np.abs(np.array(out) - np.array(ref))) == pytest.approx(
        0.0, abs=1e-10)


def test_unusable_name_still_raises():
    """Widening the pair branch must not turn it into an accept-anything:
    a genuinely unusable name= still raises, and still says what the two
    accepted forms are."""
    sc = _heisenberg()
    for bad in [(1.0, 2.0), ("Sz", "Sz"), 3.0, (sc.Sz[0],)]:
        with pytest.raises(ValueError) as excinfo:
            sc.get_dynamical_correlator(name=bad, es=ES, delta=DELTA,
                                        mode="ED", submode="ED")
        assert "documented correlator string" in str(excinfo.value)


def test_mixed_pair_is_accepted():
    """A MultiOperator on one side and a compiled operator on the other:
    accepted per-element, which is what the ED path did before name=
    resolution was centralized."""
    sc = _heisenberg()
    ref = _correlator(sc, sc.Sz[0], sc.Sz[1], mode="ED", submode="ED")
    sc = _heisenberg()
    out = _correlator(sc, sc.toMPO(sc.Sz[0], mode="ED"), sc.Sz[1],
                      mode="ED", submode="ED")
    assert np.max(np.abs(np.array(out) - np.array(ref))) == pytest.approx(
        0.0, abs=1e-10)
