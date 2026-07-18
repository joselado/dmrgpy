# Self-check for dmrgpy.pyitensor Phase 1 (Index/ITensor core). Not a
# physics example like the rest of examples/ -- there's no model here yet,
# just the tensor engine itself -- but it follows the same "run this script
# directly to check the change works" convention described in CLAUDE.md.
import numpy as np

from dmrgpy.pyitensor import (Index, ITensor, prime, noPrime, mapPrime,
                               swapPrime, dag, commonIndex, delta, svd)


def check(name, cond):
    status = "OK" if cond else "FAIL"
    print("[{}] {}".format(status, name))
    if not cond:
        raise SystemExit("selftest failed: " + name)


def test_contraction():
    i = Index(3, "Site")
    j = Index(4, "Link")
    k = Index(5, "Link")
    A = ITensor((i, j), np.random.rand(3, 4) + 1j * np.random.rand(3, 4))
    B = ITensor((j, k), np.random.rand(4, 5) + 1j * np.random.rand(4, 5))
    C = A * B
    check("contraction result indices", set(C.inds) == {i, k})
    expect = A.transpose_to((i, j)) @ B.transpose_to((j, k))
    got = C.transpose_to((i, k))
    check("contraction matches manual matmul", np.allclose(got, expect))


def test_scalar_contraction():
    i = Index(3, "Site")
    A = ITensor((i,), np.array([1, 2j, 3], dtype=complex))
    B = ITensor((i,), np.array([1, 1, 1], dtype=complex))
    s = (A * B).scalar()
    check("full contraction to scalar", np.isclose(s, 1 + 2j + 3))


def test_addition():
    i = Index(2, "Site")
    j = Index(2, "Link")
    A = ITensor((i, j), np.arange(4).reshape(2, 2).astype(complex))
    B = ITensor((j, i), np.arange(4).reshape(2, 2).astype(complex))
    C = A + B
    check("addition with transposed operand", np.allclose(
        C.transpose_to((i, j)), A.array + B.array.T))
    k = Index(2, "Link")
    D = ITensor((i, k), np.zeros((2, 2), dtype=complex))
    try:
        A + D
        ok = False
    except ValueError:
        ok = True
    check("addition rejects mismatched index structure", ok)


def test_prime_and_dag():
    i = Index(2, "Site")
    j = Index(3, "Link")
    A = ITensor((i, j), np.random.rand(2, 3) + 1j * np.random.rand(2, 3))

    Ap = prime(A, "Site")
    check("prime by tag only bumps tagged index", Ap.inds[0].plev == 1 and Ap.inds[1].plev == 0)
    check("prime preserves identity (same id)", Ap.inds[0].id == i.id)

    # Priming only "Site" leaves the "Link" leg (j) untouched, so A*Ap still
    # contracts over j (same id+plev on both sides) but can no longer
    # contract over i: A's i is at plev 0, Ap's i is at plev 1, so they're
    # distinguishable and both survive as free legs of the result.
    outer = A * Ap
    check("priming a leg makes it distinguishable from its own unprimed self",
          outer.rank == 2 and set(outer.inds) == {i, i.prime(1)})

    back = noPrime(Ap, "Site")
    check("noPrime undoes prime", set(back.inds) == {i, j})

    B = mapPrime(prime(A, "Site"), 1, 2, "Site")
    check("mapPrime moves 1->2", B.inds[0].plev == 2)

    C = swapPrime(prime(A, "Site"), 0, 1)
    check("swapPrime(0,1) on an already-primed tensor: Link(0)->1, Site(1)->0",
          C.inds[0].plev == 0 and C.inds[1].plev == 1)

    Ad = dag(A)
    check("dag conjugates values", np.allclose(Ad.array, np.conj(A.array)))
    check("dag leaves indices alone", Ad.inds == A.inds)

    Aexplicit = prime(A, i)
    check("prime by explicit Index only bumps that one", Aexplicit.inds[0].plev == 1
          and Aexplicit.inds[1].plev == 0)


def test_common_index():
    i, j, k = Index(2, "Site"), Index(3, "Link"), Index(4, "Site")
    A = ITensor((i, j), np.zeros((2, 3), dtype=complex))
    B = ITensor((j, k), np.zeros((3, 4), dtype=complex))
    check("commonIndex finds the shared bond", commonIndex(A, B) == j)
    C = ITensor((k,), np.zeros((4,), dtype=complex))
    check("commonIndex returns None when there isn't one", commonIndex(A, C) is None)


def test_delta():
    i = Index(3, "Link")
    j = i.sim()
    D = delta(i, j)
    check("delta is the identity matrix", np.allclose(D.array, np.eye(3)))


def test_svd_reconstruction():
    i, j, k = Index(3, "Site"), Index(4, "Site"), Index(5, "Site")
    data = np.random.rand(3, 4, 5) + 1j * np.random.rand(3, 4, 5)
    T = ITensor((i, j, k), data)

    U, S, V, spec = svd(T, [i, j], cutoff=0.0, maxdim=None)
    recon = (U * S) * V
    check("untruncated svd reconstructs T exactly", np.allclose(
        recon.transpose_to((i, j, k)), data, atol=1e-10))
    check("spectrum accounts for all weight when untruncated",
          np.isclose(sum(spec.eigs()) + spec.truncerr, 1.0))

    U2, S2, V2, spec2 = svd(T, [i, j], cutoff=0.0, maxdim=2)
    check("maxdim actually caps the bond dimension", U2.inds[-1].dim == 2)
    check("truncation discards some weight", spec2.truncerr > 0)
    recon2 = (U2 * S2) * V2
    err = np.linalg.norm(recon2.transpose_to((i, j, k)) - data) / np.linalg.norm(data)
    check("truncated reconstruction error is small but nonzero", 0 < err < 0.5)


if __name__ == "__main__":
    test_contraction()
    test_scalar_contraction()
    test_addition()
    test_prime_and_dag()
    test_common_index()
    test_delta()
    test_svd_reconstruction()
    print("All pyitensor Phase 1 checks passed.")
