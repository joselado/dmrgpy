"""Several spins of a chain on a metal, each with its own Wilson chain
(ysr.spin_chain_on_metal with a list of spins), and a chain whose
contacted spin is an Anderson level (ysr.anderson_spin_chain): the site
orders against each other, the decoupled limits, a spectator bath that
must not move the spectrum, the one-level and Schrieffer-Wolff limits,
and TASK against ED."""
import numpy as np
import pytest

from dmrgpy import wilsonchain, ysr

_ES, _DELTA = np.linspace(-1.2, 1.2, 13), 0.06
_ED = dict(mode="ED", submode="CVM", delta=_DELTA)
_EXPAND = dict(expand_bonds=8, expand_rounds=3)


def _free_wilson_energy(n):
    eps, t = wilsonchain.wilson_chain(n)
    e = np.linalg.eigvalsh(np.diag(eps) + np.diag(t, 1) + np.diag(t, -1))
    return 2.0*np.sum(e[e < 0])


def _level_spectrum(fc, site):
    """The spectral function of the orbital on site, summed over spin,
    as a pair of ED correlators (the reference for an Anderson tip)."""
    out = 0
    for cd, c in ((fc.Cdagup[site], fc.Cup[site]),
                  (fc.Cdagdn[site], fc.Cdn[site])):
        _, x = fc.get_dynamical_correlator(name=(c, cd), es=_ES, **_ED)
        _, y = fc.get_dynamical_correlator(name=(cd, c), es=-_ES[::-1], **_ED)
        out = out + np.real(x) + np.real(y)[::-1]
    return out


def test_one_spin_in_a_list_is_the_chain_of_one_spin():
    a = ysr.spin_chain_on_metal(3, 3, site=[1], order="linear", JK=1.0)
    b = ysr.spin_chain_on_metal(3, 3, site=1, JK=1.0)
    assert (a.hamiltonian - b.hamiltonian).is_zero()
    assert a.wilson == b.wilson and a.tip_spin == b.tip_spin


def test_layouts():
    c = ysr.spin_chain_on_metal(3, 2, site=[2, 0], order="comb")
    assert c.spin_sites == [0, 3, 4] and c.contacts == [2, 0]
    assert c.impurities == [4, 0] and c.bath_sites == [[5, 6], [1, 2]]
    assert c.tip_spin == 4 and not hasattr(c, "wilson")
    lin = ysr.spin_chain_on_metal(3, 2, site=[2, 0], order="linear")
    assert lin.spin_sites == [0, 1, 2] and lin.bath_sites == [[3, 4], [5, 6]]
    assert ysr.spin_chain_on_metal(3, 2, site=[0, 2]).order == "sides"
    s = ysr.spin_chain_on_metal(4, 2, site=[0, 1, 3])
    assert s.spin_sites == [4, 5, 6, 7] # spins 0, 1 on the left, 3 on the right
    assert s.bath_sites == [[3, 2], [1, 0], [8, 9]]


def test_comb_with_an_uncoupled_spin_warns():
    """At JK = 0 a comb solve can stop far above the ground state (0.36
    on four chains), so asking for it warns; sides does not."""
    with pytest.warns(RuntimeWarning, match="zero coupling"):
        ysr.spin_chain_on_metal(3, 2, site=[0, 2], JK=[0.5, 0.0],
                                order="comb")
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ysr.spin_chain_on_metal(3, 2, site=[0, 2], JK=[0.5, 0.0])
        ysr.spin_chain_on_metal(3, 2, site=[0, 2], JK=0.5, order="comb")


def _two_spins(order, site=(0, 1), JK=1.0, version="python"):
    fc = ysr.spin_chain_on_metal(2, 3, site=list(site), Js=0.5, JK=JK,
                                 order=order, itensor_version=version)
    fc.maxm, fc.nsweeps = 256, 20
    return fc


@pytest.mark.parametrize("order,site", [("sides", (0, 1)), ("comb", (0, 1)),
                                        ("linear", (1, 0))])
def test_plain_task_on_several_baths_raises(order, site):
    """Plain TASK misses most of X^dag|GS> on several Kondo baths unless
    the tip spin is the first of site= in the linear order (Fable's
    decision of 2026-11-18): the tip on spin 0 raises in sides, comb and
    in linear with spin 0 second, before any solve."""
    fc = _two_spins(order, site)
    with pytest.raises(ValueError, match="tip spin is the first of site="):
        ysr.get_stm_spectrum(fc, _ES, delta=_DELTA, imp=fc.spin_sites[0])


def test_the_other_routes_run_on_several_baths():
    """With bond expansion, pad_bonds="full" or mode="ED" the call runs in
    any order; the expansion's missing weight is checked against
    ysr.MISSING_RAISE, which 8/3 is far below here."""
    import warnings
    ref = _two_spins("sides")
    e = ysr.get_stm_spectrum(ref, _ES, **_ED)[1]
    fc = _two_spins("sides")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = ysr.get_stm_spectrum(fc, _ES, delta=_DELTA, **_EXPAND)[1]
        b = ysr.get_stm_spectrum(fc, _ES, delta=_DELTA, pad_bonds="full")[1]
    assert np.max(np.abs(a - e)) < 1e-8 and np.max(np.abs(b - e)) < 1e-8


def test_linear_with_the_tip_first_runs_and_warns():
    """The tip spin first in the linear order: the vertex is in the
    tangent space, 1.7e-4 of the peak against ED here (the Krylov space
    leaves it, not the start), with a warning until measured further."""
    fc = _two_spins("linear", version=3)
    ref = _two_spins("linear")
    e = ysr.get_stm_spectrum(ref, _ES, **_ED)[1]
    with pytest.warns(RuntimeWarning, match="tip spin first"):
        d = ysr.get_stm_spectrum(fc, _ES, delta=_DELTA)[1]
    assert np.max(np.abs(d - e)) < 1e-3*np.max(e)


def test_a_reordered_list_is_the_same_model():
    a = ysr.spin_chain_on_metal(3, 2, site=[2, 0], JK=[0.7, 0.4], Js=0.5,
                                order="linear")
    b = ysr.spin_chain_on_metal(3, 2, site=[0, 2], JK=[0.4, 0.7], Js=0.5)
    assert a.gs_energy(mode="ED") == pytest.approx(b.gs_energy(mode="ED"),
                                                   abs=1e-12)


def test_the_missing_weight_check_raises(monkeypatch):
    """A TASK space that lost more than ysr.MISSING_RAISE of the sum rule
    raises on several baths (forced here through get_poles)."""
    from dmrgpy import taskdmrg
    real = taskdmrg.get_poles
    def lossy(*a, **k):
        p = real(*a, **k)
        p["missing_weight"] = 10*abs(complex(np.sum(p["weight"])))
        return p
    monkeypatch.setattr(taskdmrg, "get_poles", lossy)
    fc = ysr.spin_chain_on_metal(2, 2, site=[0, 1], itensor_version="python")
    fc.maxm, fc.nsweeps = 32, 8
    with pytest.raises(ValueError, match="misses .* of\s+the sum rule"):
        ysr.get_stm_spectrum(fc, _ES, delta=_DELTA, nkry=10, expand_bonds=2)


def test_one_bath_and_levels_do_not_raise():
    import warnings
    es = np.linspace(-1.0, 1.0, 5)
    one = ysr.spin_chain_on_metal(2, 2, site=1, itensor_version="python")
    lev = ysr.anderson_spin_chain(2, 2, site=[0, 1], itensor_version="python")
    for fc in (one, lev):
        fc.maxm, fc.nsweeps = 16, 6
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ysr.get_stm_spectrum(one, es, delta=0.1, nkry=10)
    with pytest.warns(RuntimeWarning, match="levels were not measured"):
        ysr.get_stm_spectrum(lev, es, delta=0.1, nkry=10)


def test_imp_must_be_a_spin():
    fc = ysr.spin_chain_on_metal(3, 2, site=[0], order="sides")
    assert fc.spin_sites == [2, 3, 4]
    with pytest.raises(ValueError, match="is not a spin of this chain"):
        ysr.tunneling_operators(fc, imp=0)
    with pytest.raises(ValueError, match="several sites at once"):
        two = ysr.spin_chain_on_metal(3, 2, site=[0, 2])
        ysr.tunneling_operators(two, imp=two.impurities)


def test_q_is_refused_on_a_kondo_spin_by_every_route():
    fc = ysr.spin_chain_on_metal(1, 4, JK=0.4, Lambda=4.0)
    with pytest.raises(ValueError, match="q= is"):
        ysr.get_stm_spectrum(fc, _ES, submode="NRG", head="exact", q=0.3)


def test_the_two_orders_agree_on_ed():
    """Energy and the dI/dV under each spin on the metal, comb against
    linear: 4e-16 and 4e-14 of a peak of 0.44."""
    out = {}
    for order in ysr.METAL_ORDERS:
        fc = ysr.spin_chain_on_metal(3, 2, site=[0, 2], order=order,
                                     Js=0.5, JK=1.0)
        e = fc.gs_energy(mode="ED")
        d = [ysr.get_stm_spectrum(fc, _ES, imp=p, **_ED)[1]
             for p in fc.impurities]
        out[order] = (e, d)
    (ea, da), (eb, db) = out["comb"], out["linear"]
    assert ea == pytest.approx(eb, abs=1e-12)
    for x, y in zip(da, db):
        assert np.max(np.abs(x - y)) < 1e-10
    assert np.max(np.abs(da[0] - da[1])) < 1e-10 # the chain is a mirror


@pytest.mark.parametrize("order", ysr.METAL_ORDERS)
def test_decoupled_energy_is_heisenberg_plus_two_free_metals(order):
    fc = ysr.spin_chain_on_metal(2, 3, site=[0, 1], order=order, Js=0.5,
                                 JK=0.0)
    exact = -0.75*0.5 + 2*_free_wilson_energy(3)
    assert fc.gs_energy(mode="ED") == pytest.approx(exact, abs=1e-10)


def _spectator_pair(version, order):
    fc = ysr.spin_chain_on_metal(2, 3, site=[0, 1], order=order, Js=0.0,
                                 JK=0.8, itensor_version=version)
    one = ysr.spin_chain_on_metal(1, 3, JK=0.8)
    return fc, ysr.get_stm_spectrum(one, _ES, **_ED)[1]


@pytest.mark.parametrize("order", ysr.METAL_ORDERS)
def test_a_spectator_bath_does_not_move_the_spectrum_ed(order):
    """Two spins with Js=0, each on its own metal: the dI/dV under one of
    them is that of one spin on one metal, exactly."""
    fc, ref = _spectator_pair(3, order)
    for p in fc.impurities:
        d = ysr.get_stm_spectrum(fc, _ES, imp=p, **_ED)[1]
        assert np.max(np.abs(d - ref)) < 1e-10


@pytest.mark.parametrize("version", ["python", 3])
def test_a_spectator_bath_does_not_move_the_spectrum_task(version):
    fc, ref = _spectator_pair(version, "comb")
    fc.maxm, fc.nsweeps = 64, 20
    for p in fc.impurities:
        d = ysr.get_stm_spectrum(fc, _ES, imp=p, delta=_DELTA, **_EXPAND)[1]
        assert np.max(np.abs(d - ref)) < 1e-6*np.max(ref)


@pytest.mark.parametrize("version,order", [("python", "sides"), (3, "sides"),
                                           (3, "comb"), (3, "linear")])
def test_two_baths_task_against_ed(version, order):
    """Two spins, each on its own metal of three Wilson sites: 8 sites,
    saturated at maxm=256. Without bond expansion the one-site tangent
    space misses part of X^dag|GS> even so (each bath keeps its own
    charge, which caps the Schmidt rank at the bonds the vertex needs),
    up to 9.6e-2 of the peak; with expand_bonds=8, expand_rounds=3 sides
    and comb are exact under both tips. Linear keeps 1.7e-4 of the peak
    under the first spin, the same at nkry=400 and with expand_bonds=16
    (pad_bonds="full" is exact, 3e-12): the expansion is not a complete
    cure in every order."""
    kw = dict(site=[0, 1], Js=0.5, JK=1.0, order=order)
    ref = ysr.spin_chain_on_metal(2, 3, **kw)
    fc = ysr.spin_chain_on_metal(2, 3, itensor_version=version, **kw)
    fc.maxm, fc.nsweeps = 256, 20
    for p in fc.impurities:
        d = ysr.get_stm_spectrum(fc, _ES, imp=p, delta=_DELTA, **_EXPAND)[1]
        e = ysr.get_stm_spectrum(ref, _ES, imp=p, **_ED)[1]
        tol = 3e-4*np.max(e) if order == "linear" else 1e-8
        assert np.max(np.abs(d - e)) < tol


def test_tip_with_no_metal_under_it_raises():
    fc = ysr.spin_chain_on_metal(3, 2, site=[0, 2])
    with pytest.raises(ValueError, match="no Wilson chain of its own"):
        ysr.tunneling_operators(fc, imp=fc.spin_sites[1])
    # by hand it is allowed, as on a chain with one metal
    assert len(ysr.tunneling_operators(fc, imp=fc.spin_sites[1],
                                       site=fc.bath_sites[0][0])) == 2


def test_the_default_tip_is_the_first_spin_of_the_list():
    fc = ysr.spin_chain_on_metal(3, 2, site=[2, 0], Js=0.4, JK=[1.0, 0.5])
    a = ysr.get_stm_spectrum(fc, _ES, **_ED)[1]
    b = ysr.get_stm_spectrum(fc, _ES, imp=fc.impurities[0], **_ED)[1]
    c = ysr.get_stm_spectrum(fc, _ES, imp=fc.impurities[1], **_ED)[1]
    assert np.max(np.abs(a - b)) < 1e-12
    assert np.max(np.abs(a - c)) > 1e-2*np.max(a)


def test_wrong_arguments_raise():
    with pytest.raises(ValueError, match="one value per spin"):
        ysr.spin_chain_on_metal(3, 2, site=[0, 2], JK=[1.0, 2.0, 3.0])
    with pytest.raises(ValueError, match="twice"):
        ysr.spin_chain_on_metal(3, 2, site=[0, 0])
    with pytest.raises(ValueError, match="site must be one of the spins"):
        ysr.spin_chain_on_metal(3, 2, site=[0, 3])
    with pytest.raises(ValueError, match="order must be one of"):
        ysr.spin_chain_on_metal(3, 2, site=[0, 2], order="fold")
    with pytest.raises(ValueError, match="order= places several"):
        ysr.spin_chain_on_metal(3, 2, site=1, order="comb")
    fc = ysr.spin_chain_on_metal(3, 2, site=[0, 2])
    with pytest.raises(ValueError, match="runs one\\s+Wilson chain"):
        ysr.get_stm_spectrum(fc, _ES, submode="NRG")
    with pytest.raises(ValueError, match="q= is the direct tunneling"):
        ysr.tunneling_operators(fc, q=0.3)


def test_one_anderson_level_is_anderson_chain():
    a = ysr.anderson_spin_chain(1, 3, ed=-0.3, U=1.0, tD=0.3)
    b = wilsonchain.anderson_chain(3, ed=-0.3, U=1.0, tD=0.3)
    assert a.gs_energy(mode="ED") == pytest.approx(b.gs_energy(mode="ED"),
                                                   abs=1e-12)
    d = ysr.get_stm_spectrum(a, _ES, **_ED)[1]
    assert np.max(np.abs(d - _level_spectrum(b, 0))) < 1e-10


def test_anderson_vertex():
    fc = ysr.anderson_spin_chain(2, 2, U=1.0, tD=0.3)
    p, w0 = fc.anderson[0], fc.bath_sites[0][0]
    X = ysr.tunneling_operators(fc)
    assert (X[0] - fc.Cdagup[p]).is_zero() and (X[1] - fc.Cdagdn[p]).is_zero()
    X = ysr.tunneling_operators(fc, q=0.25)
    assert (X[0] - fc.Cdagup[p] - 0.25*fc.Cdagup[w0]).is_zero()
    with pytest.raises(ValueError, match="takes q="):
        ysr.tunneling_operators(fc, U=0.2)
    assert fc.penalized == [0] and fc.J_chain is None
    assert fc.ed[0] == -0.5 # the particle-hole symmetric default


def test_anderson_level_at_large_u_is_the_kondo_spin():
    """Schrieffer-Wolff: at JK = 8 tD^2/U fixed, the low-energy levels of
    the Anderson chain approach those of spin_chain_on_metal and its dI/dV
    approaches 64 (tD/U)^2 times the Kondo one, both as 1/U. Measured
    (2 spins, 3 Wilson sites, |w| < 0.6): the ratio over 64 (tD/U)^2 is
    0.919 and 0.957 at U = 32 and 64, the shape 15% and 7.9% off, the
    levels 8.5e-3 and 4.5e-3 off."""
    JK, es = 0.6, np.linspace(-0.6, 0.6, 121)
    kw = dict(mode="ED", submode="CVM", delta=0.02)
    k = ysr.spin_chain_on_metal(2, 3, site=1, Js=0.3, JK=JK)
    Ak = ysr.get_stm_spectrum(k, es, **kw)[1]
    ek = np.linalg.eigvalsh(k.get_ED_obj().get_operator(k.hamiltonian).todense())
    low = ek[ek < ek[0] + 1.0] - ek[0]
    err = []
    for U in (32.0, 64.0):
        tD = np.sqrt(JK*U/8)
        a = ysr.anderson_spin_chain(2, 3, site=1, Js=0.3, U=U, tD=tD)
        Ad = ysr.get_stm_spectrum(a, es, **kw)[1]
        ea = np.linalg.eigvalsh(a.get_ED_obj().get_operator(
            a.hamiltonian).todense())
        r = np.trapezoid(Ad, es)/np.trapezoid(Ak, es)/(64*tD**2/U**2)
        err.append((1 - r, np.max(np.abs(Ad/(64*tD**2/U**2) - Ak))/np.max(Ak),
                    np.max(np.abs(ea[:len(low)] - ea[0] - low))))
    for x, y in zip(err[0], err[1]):
        assert 0 < y < 0.6*x # every error falls about as 1/U
    assert err[1][2] < 5e-3


@pytest.mark.parametrize("version", ["python", 3])
@pytest.mark.parametrize("ns,ed,q", [(1, -0.25, 0.0), (2, -0.15, 0.2)])
def test_anderson_head_with_a_complete_basis_is_exact(version, ns, ed, q):
    """The NRG route on an Anderson level: with pad_bonds="full" the
    hand-over basis is the whole block, and with keep=None the tail keeps
    every state, so head plus tail is the exact spectrum of the chain,
    away from particle-hole symmetry (where the charge of the ground state
    is found by head_state's Nf="auto") and with the direct path q."""
    from dmrgpy import nrgmps, nrgspectral
    eps, t = wilsonchain.wilson_chain(4, Lambda=4.0, band=(-1.0, 1.0), z=0.25)
    kw = dict(Js=0.3, ed=ed, U=0.6, tD=0.3, wilson=(eps, t), Lambda=4.0)
    fc = ysr.anderson_spin_chain(ns, 4, itensor_version=version, **kw)
    fc.maxm, fc.nsweeps = 64, 20
    c, e0 = nrgmps.head_state(fc)
    assert e0 == pytest.approx(ysr.anderson_spin_chain(ns, 4, **kw)
                               .gs_energy(mode="ED"), abs=1e-9)
    sh, ex = nrgmps.handover_shell(c, 0, pad_bonds="full", Lambda=4.0, q=q)
    r = nrgspectral.spectral_poles(sh, eps[1:], t, keep=None,
                                   bins_per_decade=None)
    d = 0.05
    A = np.array([np.sum(r.weight*d/np.pi/((w - r.omega)**2 + d**2))
                  for w in _ES])
    ref = ysr.anderson_spin_chain(ns, 4, **kw)
    B = ysr.get_stm_spectrum(ref, _ES, q=q, mode="ED", submode="CVM",
                             delta=d)[1]
    assert np.max(np.abs(A - B)) < 1e-9*np.max(B)


@pytest.mark.parametrize("version", ["python", 3])
def test_anderson_task_against_ed(version):
    """An Anderson level in the middle of three spins, off particle-hole
    symmetry, with the direct path q: TASK against ED at saturation."""
    kw = dict(Js=0.4, ed=-0.3, U=1.0, tD=0.35)
    ref = ysr.anderson_spin_chain(3, 3, site=1, **kw)
    fc = ysr.anderson_spin_chain(3, 3, site=1, itensor_version=version, **kw)
    fc.maxm, fc.nsweeps = 64, 20
    e = ysr.get_stm_spectrum(ref, _ES, q=0.2, **_ED)[1]
    d = ysr.get_stm_spectrum(fc, _ES, q=0.2, delta=_DELTA)[1]
    assert np.max(np.abs(d - e)) < 1e-8


def test_the_tip_takes_its_own_bath_on_a_surrogate_chain():
    """_default_site resolves the bath of the tip spin on any chain with
    several baths; on a surrogate chain of two spins it took the first
    spin's bath for the second spin before."""
    from dmrgpy.ysrtk import surrogate
    fc = surrogate.surrogate_kondo_chain(N=2, J=0.3, Delta=0.1, L=2)
    a = ysr.tunneling_operators(fc, imp=fc.impurities[1])
    b = ysr.tunneling_operators(fc, imp=fc.impurities[1],
                                site=fc.bath_sites[1][0])
    for x, y in zip(a, b):
        assert (x - y).is_zero()


def test_head_state_with_a_fixed_sz_scans_the_charges_of_its_parity():
    """Nf="auto" on an Anderson chain with 2Sz given solves only the
    charges of the parity of 2Sz, instead of raising on the others."""
    from dmrgpy import nrgmps
    eps, t = wilsonchain.wilson_chain(3, Lambda=4.0, band=(-1.0, 1.0))
    fc = ysr.anderson_spin_chain(1, 3, U=0.6, tD=0.3, wilson=(eps, t),
                                 itensor_version="python")
    fc.maxm, fc.nsweeps = 32, 10
    for sz in (0, 1):
        c, e = nrgmps.head_state(fc, Sz=sz)
        assert (c.head_sector["Nf"] - sz) % 2 == 0
