"""edtk.dynamics.solve_cv_pair, the correction vector at +delta and at
-delta from one conjugate gradient solve, against the two solve_cv calls
it replaced: the same numbers to the last bit, on a Hermitian matrix and
a vector that are complex, as the ED route of the dynamical correlator
hands them."""
import numpy as np
import scipy.sparse as sp

from dmrgpy.edtk import dynamics


def test_one_solve_gives_both_signs_to_the_last_bit():
    rng = np.random.default_rng(7)
    n = 60
    m = sp.random(n, n, density=0.1, random_state=3) \
        + 1j*sp.random(n, n, density=0.1, random_state=4)
    h0 = sp.csr_matrix(m + m.conj().T)
    wf0 = rng.normal(size=n) + 1j*rng.normal(size=n)
    wf0 = wf0/np.linalg.norm(wf0)
    A = sp.csr_matrix(sp.random(n, n, density=0.2, random_state=5))
    B = sp.csr_matrix(A.T)
    for w in (-1.3, 0.0, 0.7):
        for delta in (0.05, 0.3):
            o1, o2 = dynamics.solve_cv_pair(h0, wf0, A, B, w, delta=delta)
            r1 = dynamics.solve_cv(h0, wf0, A, B, w, delta=delta)
            r2 = dynamics.solve_cv(h0, wf0, A, B, w, delta=-delta)
            assert o1 == r1 and o2 == r2
