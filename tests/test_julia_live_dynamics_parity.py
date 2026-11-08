"""itensor_version="julia_live" against mode="ED" and "python" on the
dynamical-correlator dispatch, the KPM moments and distributions, real-time
evolution, clones and the short-chain fallbacks.

julia_live used to have its own dispatcher in mpsjulialive/dynamics.py,
which refused submode="TD" and "ROOTN" although both run on its MPS
algebra, swallowed unknown KPM keywords, ignored kpm_energy_truncate, had
no get_dynamical_correlator_moments and no get_distribution. It now takes
the dispatch of the session backends (dynamics.py), and each submode that
lacks julia_live support refuses on its own.

Every test needs a Julia toolchain; the module is collected without one
and every test skips.
"""
import numpy as np
import pytest

from dmrgpy import spinchain, timedependent
from _helpers import julia_available

JULIA = "julia_live"


@pytest.fixture(autouse=True)
def _require_julia():
    ok, reason = julia_available()
    if not ok:
        pytest.skip("requires a working juliacall/Julia toolchain: %s"
                    % reason)


def heisenberg(n, version=JULIA, maxm=40, nsweeps=10):
    sc = spinchain.Spin_Chain([2] * n, itensor_version=version)
    sc.maxm = maxm
    sc.nsweeps = nsweeps
    h = sum(sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
            + sc.Sz[i] * sc.Sz[i + 1] for i in range(n - 1))
    # a small field on site 0 lifts the degeneracy of the ground state
    sc.set_hamiltonian(h + 0.1 * sc.Sz[0])
    return sc


def ed_spectrum(sc, version):
    """The exact spectrum of sc's Hamiltonian, through a chain of the same
    sites on mode="ED"."""
    hm = sc.get_ED_obj().get_hamiltonian()
    hm = hm.toarray() if hasattr(hm, "toarray") else np.asarray(hm)
    return np.linalg.eigvals(hm) if version == "nh" else np.linalg.eigvalsh(hm)


# ---------------------------------------------------------------- TD

PAIRS = {"SzSz": lambda sc: (sc.Sz[0], sc.Sz[0]),
         "SxSy": lambda sc: (sc.Sx[0], sc.Sy[1])}


@pytest.mark.parametrize("method", ["TDVP", "TEBD", "TDVP_GSE"])
def test_td_matches_ed_and_python(method):
    """6 sites, delta=0.2: the TD transform against the exact Lehmann
    density and against "python" TD, at the finite-window floor."""
    es = np.linspace(-0.5, 3.0, 36)
    delta = 0.2
    for key, pair in PAIRS.items():
        out = {}
        for version in (JULIA, "python"):
            sc = heisenberg(6, version)
            sc.tevol_method = method
            _, y = sc.get_dynamical_correlator(submode="TD", name=pair(sc),
                                               es=es, delta=delta)
            out[version] = np.asarray(y)
        ed = heisenberg(6, "python")
        _, yed = ed.get_dynamical_correlator(mode="ED", submode="ED",
                                             name=pair(ed), es=es, delta=delta)
        yed = np.asarray(yed)
        peak = np.max(np.abs(yed))
        # TEBD carries its second-order Trotter error at dt=0.1 on top of
        # the window floor, about 7e-3 of the peak of the (Sx0,Sy1) line
        floor = 1e-2 if method == "TEBD" else 2e-3
        assert np.max(np.abs(out[JULIA] - yed)) < floor * peak, key
        assert np.max(np.abs(out[JULIA] - out["python"])) < 2e-3 * peak, key


# ------------------------------------------------------------- ROOTN

def test_rootn_matches_ed_rootn():
    es = np.linspace(0.0, 2.5, 11)
    sc = heisenberg(4)
    kw = dict(name=(sc.Sz[0], sc.Sz[0]), es=es, delta=0.2, N=6, nkry=16)
    _, yjl = sc.get_dynamical_correlator(submode="ROOTN", **kw)
    _, yed = sc.get_dynamical_correlator(mode="ED", submode="ROOTN", **kw)
    assert np.max(np.abs(np.asarray(yjl) - np.asarray(yed))) < 1e-6


# --------------------------------------------------------------- KPM

def test_kpm_moments_are_the_moments_of_submode_kpm(monkeypatch):
    """Bit for bit: the band edge is pinned, since a julia_live solve
    starts from a random MPS and the bound would differ between calls."""
    from dmrgpy.mpsjulialive import dynamics as jldyn
    sc = heisenberg(6)
    emax = float(np.max(ed_spectrum(sc, "h")))
    monkeypatch.setattr(jldyn, "_max_energy_bound", lambda self, H: emax)
    seen = []
    orig = jldyn.dynamical_correlator_from_moments

    def spy(mus, *args, **kwargs):
        seen.append(np.array(mus))
        return orig(mus, *args, **kwargs)
    monkeypatch.setattr(jldyn, "dynamical_correlator_from_moments", spy)
    name = (sc.Sz[0], sc.Sz[0])
    mus, emin, emax2, scale, n, delta = sc.get_dynamical_correlator_moments(
        name=name, delta=0.2)
    assert len(mus) == n
    assert emax2 == emax
    es = np.linspace(0.0, 3.0, 31)
    _, y = sc.get_dynamical_correlator(submode="KPM", name=name, es=es,
                                       delta=0.2)
    assert np.array_equal(seen[-1], mus)
    # and the curve is ED's, at the same band edges
    _, yed = sc.get_dynamical_correlator(mode="ED", submode="KPM", name=name,
                                         es=es, delta=0.2)
    assert np.max(np.abs(np.asarray(y) - np.asarray(yed))) < 1e-3 * np.max(
        np.abs(yed))


def test_kpm_typo_keyword_raises_type_error():
    sc = heisenberg(4)
    with pytest.raises(TypeError):
        sc.get_dynamical_correlator(submode="KPM", name=(sc.Sz[0], sc.Sz[0]),
                                    es=np.linspace(0, 1, 3), dleta=0.1)
    with pytest.raises(TypeError):
        sc.get_dynamical_correlator_moments(name=(sc.Sz[0], sc.Sz[0]),
                                            nsweeps=3)


def test_kpm_energy_truncate_raises_naming_the_backends():
    sc = heisenberg(4)
    sc.kpm_energy_truncate = True
    with pytest.raises(NotImplementedError, match="itensor_version=3"):
        sc.get_dynamical_correlator(submode="KPM", name=(sc.Sz[0], sc.Sz[0]),
                                    es=np.linspace(0, 1, 3), delta=0.1)


def _line_weights(x, y, vals):
    x, y = np.asarray(x), np.real(np.asarray(y))
    return [np.trapezoid(y[np.abs(x - v) < 0.5], x[np.abs(x - v) < 0.5])
            for v in vals]


@pytest.mark.parametrize("scale", [2.0, None])
def test_distribution_matches_ed(scale):
    """P(x) of X = Sz_0+Sz_1 in the ground state, 6 sites: the weight of
    each of the lines x=-1,0,1 against ED. scale=None takes the operator's
    spectral edges from solves on clones."""
    sc = heisenberg(6)
    X = sc.Sz[0] + sc.Sz[1]
    xs = np.linspace(-1.5, 1.5, 601)
    x, y = sc.get_distribution(X=X, delta=0.05, scale=scale, xs=xs)
    xe, ye = sc.get_distribution(mode="ED", X=X, delta=0.05, scale=scale,
                                 xs=xs)
    wj = _line_weights(x, y, [-1.0, 0.0, 1.0])
    we = _line_weights(xe, ye, [-1.0, 0.0, 1.0])
    assert np.allclose(wj, we, atol=5e-3)
    assert sum(wj) == pytest.approx(1.0, abs=1e-2)


def test_distribution_moments_against_python():
    """The same moments on julia_live and "python", at a fixed scale; the
    6-site ground state is nondegenerate and converged on both."""
    out, outab = [], []
    for version in (JULIA, "python"):
        sc = heisenberg(6, version)
        X = sc.Sz[0] + sc.Sz[1]
        mus, _shift, _scale = sc.get_distribution_moments(X=X, delta=0.1,
                                                          scale=2.0)
        out.append(np.asarray(mus))
        # and between A|GS> and B|GS>, applied as A*wf on julia_live
        mus, _shift, _scale = sc.get_distribution_moments(
            X=X, A=sc.Sz[2], B=sc.Sx[2] + sc.Sz[3], delta=0.1, scale=2.0)
        outab.append(np.asarray(mus))
    assert len(out[0]) == len(out[1])
    assert np.max(np.abs(out[0] - out[1])) < 1e-6
    assert np.max(np.abs(outab[0] - outab[1])) < 1e-6


# --------------------------------------------- refusals and dispatch

def test_submodes_without_julia_support_refuse_cleanly():
    """No JuliaError, no AttributeError: each submode's own check."""
    sc = heisenberg(4)
    name = (sc.Sz[0], sc.Sz[0])
    es = np.linspace(0, 1, 3)
    from dmrgpy import sectordc, taskdmrg
    if JULIA not in taskdmrg.BACKENDS:
        with pytest.raises(NotImplementedError):
            sc.get_dynamical_correlator(submode="TASK", name=name, es=es)
    try:
        sectordc._check_backend(sc)
    except NotImplementedError:
        with pytest.raises(NotImplementedError):
            sc.get_dynamical_correlator(submode="SECTOR", name=name, es=es)


def test_hermitian_only_submodes_refuse_a_non_hermitian_h():
    sc = heisenberg(4)
    sc.set_hamiltonian(sc.hamiltonian + 0.3j * sc.Sz[0])
    for submode in ("TD", "ROOTN", "TDZ"):
        with pytest.raises(NotImplementedError, match="Hermitian"):
            sc.get_dynamical_correlator(submode=submode,
                                        name=(sc.Sz[0], sc.Sz[0]),
                                        es=np.linspace(0, 1, 3))


def _applyinverse_available(sc):
    try:
        sc.applyinverse(sc.hamiltonian + 3.0, sc.get_gs())
    except (TypeError, AttributeError, NotImplementedError) as exc:
        pytest.skip("needs applyinverse on julia_live (package A): %s"
                    % str(exc)[:80])


@pytest.mark.parametrize("submode", ["CVM_explicit", "CVMimag"])
def test_explicit_resolvent_submodes_match_ed(submode):
    sc = heisenberg(4)
    _applyinverse_available(sc)
    es = np.linspace(0.2, 2.0, 7)
    name = (sc.Sz[0], sc.Sz[0])
    _, yjl = sc.get_dynamical_correlator(submode=submode, name=name, es=es,
                                         delta=0.2)
    _, yed = sc.get_dynamical_correlator(mode="ED", submode="ED", name=name,
                                         es=es, delta=0.2)
    tol = 1e-5 if submode == "CVM_explicit" else 5e-2
    assert np.max(np.abs(np.asarray(yjl).real - np.asarray(yed).real)) < (
        tol * np.max(np.abs(yed)))


# --------------------------------------------------------- real time

def test_evolve_and_measure_starts_from_the_ground_state():
    """A fresh chain, no wf=: the ground state, as with wf=get_gs(), and
    the ED evolution."""
    sc = heisenberg(4)
    sc.tevol_method = "TDVP"
    _, c0 = timedependent.evolve_and_measure(sc, operator=sc.Sz[1], nt=20,
                                             dt=0.05)
    _, c1 = timedependent.evolve_and_measure(sc, operator=sc.Sz[1], nt=20,
                                             dt=0.05, wf=sc.get_gs())
    _, ce = timedependent.evolve_and_measure(sc, operator=sc.Sz[1], nt=20,
                                             dt=0.05, mode="ED")
    assert np.max(np.abs(np.asarray(c0) - np.asarray(c1))) < 1e-6
    assert np.max(np.abs(np.asarray(c0) - np.asarray(ce))) < 1e-6


# ------------------------------------------- short chains, NH n=1, S

def test_one_site_chain_falls_back_to_ed():
    sc = spinchain.Spin_Chain([2], itensor_version=JULIA)
    sc.set_hamiltonian(0.8 * sc.Sz[0] + 0.3 * sc.Sx[0])
    assert sc.get_mode() == "ED"
    assert sc.gs_energy() == pytest.approx(sc.gs_energy(mode="ED"), abs=1e-12)
    assert sc.gs_energy() == pytest.approx(-0.5 * np.hypot(0.8, 0.3),
                                           abs=1e-12)
    with pytest.raises(RuntimeError, match="n=1"):
        sc.gs_energy_generalized(sc.Sz[0] * 0 + 1.0)


def test_non_hermitian_excited_states_n1_is_the_ground_state():
    sc = heisenberg(4)
    sc.set_hamiltonian(sc.hamiltonian + 0.3j * sc.Sz[0])
    es, wfs = sc.get_excited_states(n=1)
    e0 = sc.gs_energy()
    ev = ed_spectrum(sc, "nh")
    ref = ev[np.argmin(ev.real)]
    assert es[0] == pytest.approx(e0, abs=1e-6)
    assert es[0] == pytest.approx(ref, abs=1e-6)
    assert len(wfs) == 1


def test_bond_entropy_of_a_scaled_state():
    sc = heisenberg(8)
    wf = sc.get_gs()
    s1 = sc.get_bond_entropy(wf, 3, 4)
    s2 = sc.get_bond_entropy(2.0 * wf, 3, 4)
    assert s2 == pytest.approx(s1, abs=1e-8)
    # and the ED state's Schmidt values across the same bond
    v = np.asarray(sc.get_gs(mode="ED").v).reshape(2 ** 4, 2 ** 4)
    p = np.linalg.svd(v, compute_uv=False) ** 2
    p = p[p > 1e-14] / np.sum(p)
    assert s1 == pytest.approx(-np.sum(p * np.log(p)), abs=1e-6)


# ------------------------------------------------- clones and thermal

def test_clone_bandwidth_and_thermal_state():
    sc = heisenberg(4)
    e0 = sc.gs_energy()
    c = sc.clone()
    c.set_hamiltonian(c.hamiltonian) # a cold start: the clone solves again
    wc = c.get_gs()
    assert c.gs_energy() == pytest.approx(e0, abs=1e-8)
    assert abs(wc.dot(sc.get_gs())) == pytest.approx(1.0, abs=1e-6)
    # bandwidth() solves X and -X on clones. X has nondegenerate band
    # edges: the top of the Heisenberg chain is a ferromagnetic multiplet
    # split only by the 0.1 field, and a random-start solve of -H there
    # lands 1e-5 to 1e-4 below its top from run to run
    X = sc.hamiltonian + sum((0.5 + 0.2 * i) * sc.Sz[i] + 0.3 * sc.Sx[i]
                             for i in range(4))
    xm = np.asarray(sc.get_ED_obj().MO2matrix(X).todense())
    ev = np.linalg.eigvalsh(xm)
    assert sc.bandwidth(X) == pytest.approx(ev[-1] - ev[0], abs=1e-6)
    # the chain itself was not touched by the clone's solves
    assert sc.gs_energy() == pytest.approx(e0, abs=1e-12)
    from dmrgpy import thermal
    out = {}
    for mode in (None, "ED"):
        tc = thermal.Thermal_Spin_Chain([2] * 4, T=0.5, maxm=64,
                                        itensor_version=JULIA, mode=mode)
        h = sum(tc.Sx[i] * tc.Sx[i + 1] + tc.Sy[i] * tc.Sy[i + 1]
                + tc.Sz[i] * tc.Sz[i + 1] for i in range(3))
        tc.set_hamiltonian(h)
        wf = tc.get_gs()
        out[mode] = (wf.dot(h * wf) / wf.dot(wf)).real
    assert out[None] == pytest.approx(out["ED"], abs=1e-4)


# ---------------------------------------------- fermion Green's function

def test_get_gr_against_python():
    """fermionchain.get_gr calls kpmdmrg.get_dynamical_correlator
    directly, with name="cdc"/"ccd"; it needs a parity-odd MPO for a
    single C on julia_live."""
    from dmrgpy import fermionchain
    out = {}
    es = np.linspace(-2.0, 2.0, 41)
    for version in (JULIA, "python"):
        fc = fermionchain.Fermionic_Chain(4, itensor_version=version)
        fc.maxm = 30
        h = sum(fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
                for i in range(3))
        fc.set_hamiltonian(h + 0.1 * fc.N[0])
        if version == JULIA:
            from dmrgpy.mpsjulialive.mpo import MPO
            try:
                MPO(fc.C[1], MBO=fc) * fc.get_gs()
            except Exception as exc:
                pytest.skip("needs parity-odd MPOs on julia_live "
                            "(package B): %s" % str(exc)[:80])
        x, y = fc.get_gr(delta=0.1, es=es, i=1, j=1)
        out[version] = np.asarray(y)
    peak = np.max(np.abs(out["python"]))
    assert np.max(np.abs(out[JULIA] - out["python"])) < 5e-2 * peak
