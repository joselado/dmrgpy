"""The representative-energy (Zitko-Pruschke) Wilson mesh, mesh="zp" of
wilsonchain.star_levels, against the mean-energy mesh: the closed form,
A_Lambda = 1 on the classical YSR energy, and the band edge of the
z-averaged density of states of the BdG chain; and two riders of ysr.py,
the exact gamma=0 BCS density and the TASK warning at J=0; and the
factor of the mesh in the z-averaged spectrum of ysrtk.zaverage."""
import warnings

import numpy as np
import pytest

from dmrgpy import wilsonchain, ysr
from dmrgpy.ysrtk import zaverage


def _bdg(eps, t, Delta, v0=0.0):
    n = len(eps)
    h = np.diag(eps) + np.diag(t, 1) + np.diag(t, -1)
    H = np.block([[h, Delta*np.eye(n)], [Delta*np.eye(n), -h]])
    H[0, 0] += v0
    H[n, n] += v0
    return H


def test_mean_is_the_default_mesh():
    for L in (1.5, 2.0, 4.0):
        for z in (0.0, 0.25, 0.5):
            e0, t0 = wilsonchain.wilson_chain(20, Lambda=L, z=z)
            e1, t1 = wilsonchain.wilson_chain(20, Lambda=L, z=z, mesh="mean")
            assert np.array_equal(e0, e1) and np.array_equal(t0, t1)
        _, t = wilsonchain.wilson_chain(30, Lambda=L)
        assert np.max(np.abs(t/wilsonchain.flat_band_hoppings(30, L) - 1)) < 1e-14


def test_zp_mesh_is_checked():
    with pytest.raises(ValueError):
        wilsonchain.star_levels(mesh="ZP")
    with pytest.raises(NotImplementedError):
        wilsonchain.wilson_chain(6, dos=lambda e: 1 + e*e, mesh="zp")
    with pytest.raises(ValueError):
        ysr.kondo_chain(4, mesh="log")


@pytest.mark.parametrize("Lambda", [2.0, 4.0])
def test_zp_levels_closed_form(Lambda):
    """|dE/dz| equals the weight of each interval (eq. 35 of
    arXiv:0810.4765), the weights are those of the mean mesh, and at z=0
    every level is the mean one over A_Lambda."""
    A = ysr.lambda_factor(Lambda)
    for z in (0.1, 0.4, 0.8):
        h = 1e-6
        xp = wilsonchain.star_levels(Lambda=Lambda, z=z+h, nint=6, mesh="zp")[0]
        xm = wilsonchain.star_levels(Lambda=Lambda, z=z-h, nint=6, mesh="zp")[0]
        x, w = wilsonchain.star_levels(Lambda=Lambda, z=z, nint=6, mesh="zp")
        wm = wilsonchain.star_levels(Lambda=Lambda, z=z, nint=6)[1]
        assert np.array_equal(w, wm)
        k = np.arange(6) # the top six levels of the positive side
        dE = np.abs(xp[k] - xm[k])/(2*h)
        assert np.max(np.abs(dE/(2*w[k]) - 1)) < 1e-8 # w covers both sides
        top = (1 - Lambda**-z)/np.log(Lambda) + 1 - z
        assert x[0] == pytest.approx(top, rel=1e-14)
    x, _ = wilsonchain.star_levels(Lambda=Lambda, nint=10, mesh="zp")
    xm, _ = wilsonchain.star_levels(Lambda=Lambda, nint=10)
    assert np.max(np.abs(x[:10]*A/xm[:10] - 1)) < 1e-14
    _, t = wilsonchain.wilson_chain(20, Lambda=Lambda, mesh="zp")
    assert np.max(np.abs(t*A/wilsonchain.flat_band_hoppings(20, Lambda) - 1)) < 1e-10


@pytest.mark.parametrize("Lambda", [1.5, 2.0, 3.0])
@pytest.mark.parametrize("z", [0.0, 0.5])
def test_classical_ysr_energy_with_the_bare_field(Lambda, z):
    """On the zp mesh a static field v0 = JS/2 acts as v0 itself: the
    subgap state is at Delta (1-a^2)/(1+a^2) with a = pi rho J S/2."""
    Delta = 1e-3
    eps, t = wilsonchain.wilson_chain(80, Lambda=Lambda, z=z, mesh="zp")
    for JS in (0.3, 0.6, 1.0):
        ev = np.linalg.eigvalsh(_bdg(eps, t, Delta, v0=JS/2))
        a = np.pi*0.5*(JS/2)
        pred = Delta*abs(1 - a*a)/(1 + a*a)
        assert abs(np.min(np.abs(ev)) - pred) < 1e-3*Delta


def test_classical_ysr_energy_at_lambda_4_on_two_twists():
    """At Lambda=4 one twist misses by up to 3.3e-3 Delta, on this mesh as
    on the mean mesh at the same continuum coupling (8 twists, JS=1): the
    coarse mesh, not the chain. The two twists z=0 and 0.5 cancel it."""
    Delta = 1e-3
    for JS in (0.3, 0.6, 1.0):
        a = np.pi*0.5*(JS/2)
        pred = Delta*abs(1 - a*a)/(1 + a*a)
        e = []
        for z in (0.0, 0.5):
            eps, t = wilsonchain.wilson_chain(80, Lambda=4.0, z=z, mesh="zp")
            e.append(np.min(np.abs(np.linalg.eigvalsh(_bdg(eps, t, Delta, JS/2)))))
        assert abs(np.mean(e) - pred) < 1e-3*Delta


def _zavg(xi, mesh, nz, b, n=20, Lambda=4.0, Delta=1e-3):
    es = np.sqrt(xi**2 + Delta**2)
    A = 0.0
    for z in zaverage.twists(nz):
        eps, t = wilsonchain.wilson_chain(n, Lambda=Lambda, z=z, mesh=mesh)
        w, wt = zaverage.bdg_poles(eps, t, Delta)
        w0 = zaverage._default_omega0(w, Delta, t[-1])
        A = A + zaverage.broaden_poles(es, w, wt, b=b, omega0=w0, gap=Delta)[0]
    fac = ysr.lambda_factor(Lambda) if mesh == "mean" else 1.0
    return A/nz*fac, 0.5*ysr.bcs_dos(es, Delta)


def test_zp_density_of_states_with_no_factor():
    """The z average on the zp mesh, with no A_Lambda, is the BCS density
    of states at low energy, and carries 1/2 per side, against
    A_Lambda/2 = 0.578 for the scaled mean mesh."""
    A, ref = _zavg(np.logspace(-4, -2, 21), "zp", 4, 0.4)
    assert np.max(np.abs(A/ref - 1)) < 2e-5
    for mesh, side in (("zp", 0.5), ("mean", ysr.lambda_factor(4.0)/2)):
        c = 0.0
        for z in zaverage.twists(16):
            eps, t = wilsonchain.wilson_chain(20, Lambda=4.0, z=z, mesh=mesh)
            w, wt = zaverage.bdg_poles(eps, t, 1e-3)
            c += np.sum(wt[w > 0])/16
        fac = ysr.lambda_factor(4.0) if mesh == "mean" else 1.0
        assert c*fac == pytest.approx(side, rel=1e-10)


def test_zp_band_edge_is_flat():
    """With 64 twists and b=0.05 the zp curve is 1/2 up to 0.9 D, where
    the mean mesh piles its top levels up to 2.5 times too high."""
    xi = np.array([0.1, 0.5, 0.7, 0.9])
    A, ref = _zavg(xi, "zp", 64, 0.05)
    assert np.max(np.abs(A - ref)) < 2e-3
    A, ref = _zavg(xi, "mean", 64, 0.05)
    assert A[2] > 1.2


def test_kondo_chain_stores_the_mesh():
    fc = ysr.kondo_chain(4, J=0.3, Lambda=4.0, mesh="zp")
    assert fc.mesh == "zp" and fc.dos_factor == 1.0 and fc.J_chain == 0.3
    _, t = wilsonchain.wilson_chain(4, Lambda=4.0, mesh="zp")
    assert fc.t_last == abs(t[-1])
    fc = ysr.kondo_chain(4, J=0.3, Lambda=4.0)
    assert fc.mesh == "mean"
    assert fc.dos_factor == ysr.lambda_factor(4.0)
    _, t = wilsonchain.wilson_chain(4, Lambda=4.0)
    assert fc.t_last == abs(t[-1])
    assert ysr.kondo_chain(1, J=0.3, D=2.0).t_last == 2.0


def test_bcs_dos_at_zero_gamma_is_exact():
    """No density inside the gap, and the old form 1e-14j outside it."""
    Delta = 1e-3
    e = Delta*np.concatenate([np.linspace(-0.999, 0.999, 201),
                              1 + np.logspace(-6, 3, 200),
                              -1 - np.logspace(-6, 3, 200)])
    rho = ysr.bcs_dos(e, Delta)
    inside = np.abs(e) < Delta
    assert np.all(rho[inside] == 0.0)
    w = e[~inside] + 1e-14j
    old = np.abs(np.real(w/np.sqrt(w*w - Delta*Delta)))
    assert np.max(np.abs(rho[~inside]/old - 1)) < 1e-9
    assert np.array_equal(ysr.bcs_dos(e, -Delta), rho)
    assert np.all(ysr.bcs_dos(e, 0.0) == 1.0)
    assert ysr.bcs_dos(2*Delta, Delta) == pytest.approx(2/np.sqrt(3), rel=1e-14)
    assert ysr.bcs_dos(e, Delta, gamma=1e-5)[0] > 0.0 # Dynes unchanged


def test_task_at_zero_coupling_warns():
    """At J=0 the ground state is a product of the impurity and the bath,
    and TASK misses the spin-flip part of X|GS>; get_stm_spectrum says
    so, and the ED route, which is exact, does not warn."""
    es = np.linspace(-1.2, 1.2, 25)
    fc = ysr.kondo_chain(4, J=0.0, Delta=0.1, itensor_version="python")
    fc.maxm, fc.nsweeps = 60, 15
    with pytest.warns(RuntimeWarning, match="zaverage.bdg_chain_dos"):
        _, A = ysr.get_stm_spectrum(fc, es, delta=0.08)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        _, ed = ysr.get_stm_spectrum(fc, es, delta=0.08, mode="ED",
                                     submode="CVM")
    assert not any("J = 0 on the chain" in str(w.message) for w in rec)
    assert np.sum(A) < 0.9*np.sum(ed) # the weight TASK misses
    fc = ysr.kondo_chain(4, J=0.6, Delta=0.1, itensor_version="python")
    fc.maxm, fc.nsweeps = 60, 15
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        ysr.get_stm_spectrum(fc, es, delta=0.08)
    assert not any("J = 0 on the chain" in str(w.message) for w in rec)


@pytest.mark.parametrize("mesh", ["mean", "zp"])
def test_zaverage_spectrum_takes_the_factor_of_the_mesh(mesh):
    # the one-electron spectrum of the first bath site at J=0, through the
    # z-averaged TASK route, is the BdG density of the same mesh: the factor
    # is fc.dos_factor (A_Lambda on mean, 1 on zp), and the last hopping is
    # that of the chain's own mesh (before, zp came out A_Lambda too large)
    es = np.linspace(-1.5, 1.5, 61)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, A, info = zaverage.zaverage_stm_spectrum(
            6, es, J=0.0, Delta=0.1, nz=1, b=0.4, mesh=mesh,
            itensor_version=3, operators=lambda fc: [fc.Cdagup[1]])
    ref = zaverage.bdg_chain_dos(es, 6, 0.1, nz=1, b=0.4, mesh=mesh)
    out = np.abs(es) > 0.15
    assert info["missing_weight"][0] < 1e-10
    assert A[out] == pytest.approx(ref[out], rel=1e-8)
