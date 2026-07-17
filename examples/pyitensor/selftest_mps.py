# Self-check for dmrgpy.pyitensor Phase 4 (MPS/MPO containers + algebra:
# position/canonicalization, sum, applyMPO, nmultMPO, inner/innerC, traceC,
# randomMPS, and AutoMPO -> to_mpo). Cross-checked against dense-vector/
# matrix reference computations (testutil.py) built via an independent
# path (sequential SVD splitting, not DMRG/randomMPS), and against
# AutoMPO.dense_matrix() from Phase 3.
import numpy as np

from dmrgpy.pyitensor import (AutoMPO, applyMPO, inner, nmultMPO,
                               randomMPS, sum as mps_sum, to_mpo, traceC)
from dmrgpy.pyitensor.sites import SiteX
from testutil import dense_to_mps, mpo_to_dense, mps_to_dense


def check(name, cond):
    status = "OK" if cond else "FAIL"
    print("[{}] {}".format(status, name))
    if not cond:
        raise SystemExit("selftest failed: " + name)


def heisenberg_ampo(sites, n):
    ampo = AutoMPO(sites)
    for i in range(1, n):
        ampo.add(1.0, "Sz", i, "Sz", i + 1)
        ampo.add(0.5, "Sp", i, "Sm", i + 1)
        ampo.add(0.5, "Sm", i, "Sp", i + 1)
    return ampo


def hopping_ampo(sites, n):
    ampo = AutoMPO(sites)
    for i in range(1, n):
        ampo.add(1.0, "Cdag", i, "C", i + 1)
        ampo.add(1.0, "Cdag", i + 1, "C", i)
    ampo.add(0.4, "Cdag", 1, "C", n)
    ampo.add(0.4, "Cdag", n, "C", 1)
    return ampo


def test_random_mps_and_position():
    np.random.seed(0)
    n = 5
    sites = SiteX([2] * n)
    psi = randomMPS(sites, 8)
    check("randomMPS is normalized", np.isclose(inner(psi, psi), 1.0, atol=1e-10))
    v0 = mps_to_dense(psi)
    check("dense vector from randomMPS is normalized", np.isclose(np.vdot(v0, v0), 1.0))

    for b in [3, 1, 5, 2]:
        psi.position(b)
        check("position({}) preserves the represented state".format(b),
              np.allclose(mps_to_dense(psi), v0, atol=1e-8))


def test_sum_mps():
    np.random.seed(1)
    n = 4
    sites = SiteX([2] * n)
    a = randomMPS(sites, 4)
    b = randomMPS(sites, 4)
    da, db = mps_to_dense(a), mps_to_dense(b)
    c = mps_sum(a, b, cutoff=0.0, maxdim=1000)
    check("sum(A,B) matches dense A+B", np.allclose(mps_to_dense(c), da + db, atol=1e-8))


def test_inner():
    np.random.seed(2)
    n = 4
    sites = SiteX([0] * n)  # spinless fermion, dim 2 per site like spin-1/2
    a = randomMPS(sites, 6)
    b = randomMPS(sites, 6)
    da, db = mps_to_dense(a), mps_to_dense(b)
    check("inner(A,B) matches dense <A|B>", np.isclose(inner(a, b), np.vdot(da, db), atol=1e-8))


def test_autompo_to_mpo_spin():
    n = 4
    sites = SiteX([2] * n)
    ampo = heisenberg_ampo(sites, n)
    ref = ampo.dense_matrix()
    mpo = to_mpo(ampo, cutoff=0.0, maxdim=1000)
    check("to_mpo(Heisenberg) matches AutoMPO.dense_matrix()",
          np.allclose(mpo_to_dense(mpo), ref, atol=1e-8))


def test_autompo_to_mpo_fermion():
    n = 4
    sites = SiteX([0] * n)
    ampo = hopping_ampo(sites, n)
    ref = ampo.dense_matrix()
    mpo = to_mpo(ampo, cutoff=0.0, maxdim=1000)
    check("to_mpo(hopping, incl. long-range JW) matches AutoMPO.dense_matrix()",
          np.allclose(mpo_to_dense(mpo), ref, atol=1e-8))


def test_apply_mpo():
    np.random.seed(3)
    n = 4
    sites = SiteX([2] * n)
    ampo = heisenberg_ampo(sites, n)
    H = to_mpo(ampo, cutoff=0.0, maxdim=1000)
    psi = randomMPS(sites, 6)
    Hpsi = applyMPO(H, psi, cutoff=0.0, maxdim=1000)
    Hpsi.noPrime("Site")
    dense_ref = mpo_to_dense(H) @ mps_to_dense(psi)
    check("applyMPO(H,psi) matches dense H @ psi", np.allclose(mps_to_dense(Hpsi), dense_ref, atol=1e-8))


def test_nmult_mpo_and_trace():
    n = 3
    sites = SiteX([2] * n)
    ampo = heisenberg_ampo(sites, n)
    H = to_mpo(ampo, cutoff=0.0, maxdim=1000)
    H2 = nmultMPO(H, H.copy().prime(), cutoff=0.0, maxdim=1000)
    H2.mapPrime(2, 1)
    dense_H = mpo_to_dense(H)
    check("nmultMPO(H,prime(H)) then mapPrime(2,1) matches dense H@H",
          np.allclose(mpo_to_dense(H2), dense_H @ dense_H, atol=1e-8))

    check("traceC(H) matches dense trace", np.isclose(traceC(H), np.trace(dense_H), atol=1e-8))


def test_dense_to_mps_roundtrip():
    n = 4
    sites = SiteX([2] * n)
    v = np.random.randn(2 ** n) + 1j * np.random.randn(2 ** n)
    v /= np.linalg.norm(v)
    mps = dense_to_mps(v, sites, maxdim=1000)
    check("dense_to_mps roundtrips exactly", np.allclose(mps_to_dense(mps), v, atol=1e-8))
    check("dense_to_mps result is normalized via inner()", np.isclose(inner(mps, mps), 1.0, atol=1e-8))


if __name__ == "__main__":
    test_random_mps_and_position()
    test_sum_mps()
    test_inner()
    test_autompo_to_mpo_spin()
    test_autompo_to_mpo_fermion()
    test_apply_mpo()
    test_nmult_mpo_and_trace()
    test_dense_to_mps_roundtrip()
    print("All pyitensor Phase 4 checks passed.")
