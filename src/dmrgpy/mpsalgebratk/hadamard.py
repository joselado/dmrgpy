"""Elementwise (Hadamard) products of wavefunctions in the computational
basis, and what is built on them.

The product (psi*phi)(s) = psi(s) phi(s) is taken in the basis the site
indices carry (Sz for spins, occupations for fermions and bosons), so it is
a statement about that basis, not about the states alone. For fermions a
product of two amplitudes also depends on the sign convention of the basis
states, while |psi(s)|^2 does not.

On the MPS backends the product is computed by recursive sketched
interpolation (pyitensor/rsi.py, arXiv:2602.17974), or exactly by the
Kronecker product of the site tensors followed by an SVD compression. The
kernel runs on the site arrays, so it is the same code for
itensor_version=3 and "python"; the arrays come out of and go back into
the session (mps_site_arrays, mps_from_site_arrays). mode="ED" multiplies
the dense vectors.
"""

import numpy as np

from .. import mps


def _check_session(self, wfs, where):
    from ..mpsalgebra import _is_julia
    if any(_is_julia(wf) for wf in wfs):
        raise NotImplementedError(where+" is not available on "
                "itensor_version='julia_live'; use itensor_version=3, "
                "'python', or mode='ED'")
    session = getattr(self, "_session", None)
    if session is None or not hasattr(session, "mps_from_site_arrays"):
        raise NotImplementedError(where+" needs itensor_version=3 or "
                "'python' (got itensor_version=%r); use mode='ED' for a "
                "small chain" % (self.itensor_version,))
    if getattr(self, "conserved_sector", None):
        raise ValueError("conserved-sector mode (%s): %s works on dense "
                "sites only; clear the sector with set_conserved_sector()"
                % (self.conserved_sector, where))
    return session


def _conj_flags(conj, n):
    if isinstance(conj, (bool, np.bool_)) or conj is None:
        return [bool(conj)]*n
    flags = [bool(c) for c in conj]
    if len(flags) != n:
        raise ValueError("conj= needs one flag per wavefunction (got %d "
                         "flags for %d wavefunctions)" % (len(flags), n))
    return flags


def _arrays(session, wf, conj):
    arrs = [np.asarray(a) for a in session.mps_site_arrays(wf.cpp_handle)]
    return [a.conj() for a in arrs] if conj else arrs


def _to_mps(self, session, arrays):
    from ..pyitensor import backend as bk
    host = [np.ascontiguousarray(bk.to_host(a), dtype=complex)
            for a in arrays]
    handle = session.mps_from_site_arrays(host)
    return mps.MPS(self, cpp_handle=handle).copy()


def _compress(self, wf, maxm):
    """An exact (in bond dimension) state compressed by SVD to the chain's
    cutoff and at most maxm states per bond: the sum of the state with
    zero, which every session truncates the way sum_mps does."""
    session = self._session
    session.set_sweep_params(int(maxm), self.nsweeps, self.cutoff, self.noise)
    try:
        zero = session.scale_mps(wf.cpp_handle, 0.0)
        handle = session.sum_mps(wf.cpp_handle, zero)
    finally:
        session.set_sweep_params(self.maxm, self.nsweeps, self.cutoff,
                                 self.noise)
    return mps.MPS(self, cpp_handle=handle).copy()


def hadamard_product(self, *wfs, conj=False, method="rsi", maxm=None,
                     tol=1e-12, oversampling=8, seed=0, mode=None):
    """The elementwise product psi_1(s) psi_2(s) ... of the given states in
    the computational basis, as a state of the same chain (not normalized).

    conj: a flag, or one flag per state, conjugating that state's
        amplitudes first; hadamard_product(wf, wf, conj=(True, False)) is
        the probability |psi(s)|^2 as a state.
    method: "rsi", recursive sketched interpolation (arXiv:2602.17974),
        O(chi^3); or "direct", the exact Kronecker product compressed by
        SVD, O(chi^4) or worse, the reference for "rsi".
    maxm: the largest bond dimension of the result (default: the chain's
        maxm).
    tol: interpolation tolerance of "rsi", relative to the largest entry.
    oversampling: extra sketch columns of "rsi" above maxm/d.
    seed: seed of the random sketch of "rsi".
    mode: "ED" multiplies the dense vectors of ED states."""
    if len(wfs) < 2:
        raise ValueError("hadamard_product needs at least two wavefunctions")
    from ..mpsalgebra import wavefunction_mode
    modes = set(wavefunction_mode(wf, mode=mode) for wf in wfs)
    if len(modes) != 1:
        raise ValueError("hadamard_product: the wavefunctions belong to "
                         "different backends")
    flags = _conj_flags(conj, len(wfs))
    if modes == {"ED"}:
        from ..edtk.edchain import State
        v = None
        for wf, c in zip(wfs, flags):
            w = np.conjugate(wf.v) if c else wf.v
            v = w if v is None else v*w
        return State(v, wfs[0].MBO)
    session = _check_session(self, wfs, "hadamard_product")
    if maxm is None: maxm = self.maxm
    tts = [_arrays(session, wf, c) for wf, c in zip(wfs, flags)]
    from ..pyitensor import rsi
    if method == "rsi":
        return _to_mps(self, session, rsi.hadamard_arrays(tts, maxm,
                tol=tol, oversampling=oversampling, seed=seed))
    elif method == "direct":
        exact = _to_mps(self, session, rsi.kronecker_arrays(tts))
        return _compress(self, exact, maxm)
    raise ValueError("hadamard_product: method must be 'rsi' or 'direct' "
                     "(got %r)" % (method,))


def elementwise_map(self, wf, func, maxm=None, tol=1e-12, oversampling=8,
                    seed=0, mode=None):
    """The state with amplitudes func(psi(s)) in the computational basis,
    by the nonlinear-map extension of recursive sketched interpolation
    (arXiv:2602.17974, Sec. 6): func is applied to the sketched
    amplitudes where the product would be taken. func acts elementwise on
    a NumPy array. Exact once the bond dimension holds the result and the
    tail is short; a heuristic otherwise, as the paper says. mode="ED"
    applies func to the dense vector."""
    from ..mpsalgebra import wavefunction_mode
    if wavefunction_mode(wf, mode=mode) == "ED":
        from ..edtk.edchain import State
        return State(np.asarray(func(wf.v), dtype=complex), wf.MBO)
    session = _check_session(self, [wf], "elementwise_map")
    if maxm is None: maxm = self.maxm
    from ..pyitensor import rsi
    return _to_mps(self, session, rsi.hadamard_arrays(
            [_arrays(session, wf, False)], maxm, tol=tol,
            oversampling=oversampling, seed=seed, func=func))


def participation_entropy(self, wf=None, q=2, mode="DMRG", **kwargs):
    """The Renyi participation entropy of a state in the computational
    basis, S_q = log(sum_s p(s)^q)/(1-q) with p(s) = |psi(s)|^2/<psi|psi>,
    for an integer q >= 2. On the MPS backends p is built by
    hadamard_product (keywords as there), and sum_s p^q is the overlap of
    p with p^(q-1), the product of 2(q-1) copies of the state."""
    if int(q) != q or q < 2:
        raise ValueError("participation_entropy: q must be an integer >= 2 "
                         "(got %r)" % (q,))
    q = int(q)
    if wf is None: wf = self.get_gs(mode=mode)
    from ..mpsalgebra import wavefunction_mode
    if wavefunction_mode(wf) == "ED":
        p = np.abs(wf.v)**2
        p = p/np.sum(p)
        return float(np.log(np.sum(p**q))/(1-q))
    norm2 = wf.dot(wf).real
    p = hadamard_product(self, wf, wf, conj=(True, False), **kwargs)
    if q == 2:
        s = p.dot(p).real
    else:
        copies = [wf]*(q-1) + [wf]*(q-1)
        flags = [True]*(q-1) + [False]*(q-1)
        pq1 = hadamard_product(self, *copies, conj=flags, **kwargs)
        s = p.dot(pq1).real
    return float(np.log(s/norm2**q)/(1-q))
