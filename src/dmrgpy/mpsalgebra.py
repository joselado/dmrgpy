from . import mps
import numpy as np


def wavefunction_mode(wf,mode=None):
    """Resolve the solver that answers an MPS/operator-algebra primitive.

    For these primitives the wavefunction *is* the backend -- an
    mps.MPS lives in a DMRG session, an edtk.edchain.State is a dense
    ED vector -- so the type is what decides, not self.mode. A mode=
    given by the caller is still honoured, as a *check*: the user guide
    advertises mode= across this whole family, and silently running the
    other backend (which is what dropping it amounted to) is worse than
    saying the two disagree.

    The type test used to be `type(wf)==np.ndarray` in applyoperator/
    summps, which no ED route has produced since EDchain.get_gs started
    returning a State: both ED branches were dead and every ED call fell
    into a bare `raise`, i.e. "RuntimeError: No active exception to
    reraise". applyinverse, in the same file, already tested State."""
    from .edtk.edchain import State
    if isinstance(wf,mps.MPS): actual = "DMRG"
    elif isinstance(wf,State): actual = "ED"
    else:
        raise TypeError("unsupported wavefunction type "
                +type(wf).__name__+": expected an mps.MPS (DMRG backends) "
                "or an edtk.edchain.State (ED backend)")
    if mode is not None and mode!=actual:
        raise TypeError("mode="+repr(mode)+" was requested, but this "
                "wavefunction is a "+type(wf).__name__+", which only the "
                +actual+" backend can consume. Rebuild the wavefunction "
                "with that backend (e.g. get_gs(mode="+repr(mode)+")).")
    return actual


def exponential(self,h,wf,mode=None,**kwargs):
    """Compute exp(h)|wf>"""
    mode = wavefunction_mode(wf,mode=mode) # solver, see above
    if mode=="DMRG":
        # Gate on the chain's *numerical* Hermiticity probe, not on
        # MultiOperator.is_hermitian(). The symbolic test compares
        # h-h.get_dagger() against 0 after simplify(), which does not
        # know that get_dagger()'s factor-order reversal is a no-op for
        # factors living on different sites -- so it reports False for
        # Sx[i]*Sx[j]+Sy[i]*Sy[j]+Sz[i]*Sz[j], the single most common
        # Hamiltonian shape in this library (the same false negative
        # infinitechain.py:_check_reach_one documents for its own use).
        # Both branches then failed and control fell into an
        # uncontrolled 2-term Taylor truncation with no step
        # subdivision: 4% wrong at z=1 on a 4-site Heisenberg chain and
        # unbounded in z.
        if self.is_hermitian(h):
            return exponential_dmrg(self,h,wf,dt=1.0,**kwargs)
        elif self.is_hermitian(1j*h): # i.e. h is anti-Hermitian
            # exp(h) = exp(1j*(-1j*h)) with -1j*h Hermitian, which is
            # what exponential_dmrg's own Hermiticity check needs
            return exponential_dmrg(self,-1j*h,wf,dt=1j,**kwargs)
        else:
            raise NotImplementedError(
                "exponential() needs a Hermitian or anti-Hermitian "
                "operator on the DMRG backends (the in-process "
                "extension's custom_exp is a truncated Taylor expansion, "
                "convergent only when the step count can be set from the "
                "operator's bandwidth). This one is neither; use "
                "mode=\"ED\" on a small chain instead.")
    elif mode=="ED":
        return self.get_ED_obj().exponential(h,wf,**kwargs)
    else: raise ValueError("Unrecognized mode "+repr(mode))


def exponential_dmrg(self,h,wfa,dt=1.0,nt=1000,nt0=None):
    """Compute exp(dt*h)|wfa> via the in-process pybind11 extension
    (mpscpp2/chain_session.h's Chain::exponential_apply, a custom
    2nd-order Taylor expansion applied over nt0 sub-steps)."""
    if not self.is_hermitian(h):
        raise ValueError("exponential_dmrg needs a Hermitian operator "
                "(the sub-step count is set from its bandwidth, which is "
                "only meaningful for a real spectrum); pass -1j*h for an "
                "anti-Hermitian h, as exponential() does")
    if nt0 is None:
        nt0 = int(h.get_bandwidth(self)*nt)
        # get_bandwidth() runs its own DMRG ground-state search (see
        # bandwidth() in manybodychain.py) and can occasionally
        # underestimate a highly degenerate operator's spectral width
        # depending on the random initial wavefunction; nt0<=0 would divide
        # by a zero/negative effective dt below, so fall back to nt steps
        # (equivalent to a bandwidth of 1) rather than feeding a degenerate
        # step count into the extension.
        if nt0<1: nt0 = nt
    if not self.tevol_custom_exp:
        raise NotImplementedError(
                "tevol_custom_exp=False selects ITensor's toExpH variant, "
                "which only ever existed in the removed file-based backend; "
                "the in-process extension only implements the custom_exp "
                "(2nd-order Taylor) variant, so leave tevol_custom_exp=True")
    # Chain::exponential_apply computes exp(tau*h)|wfa>, so tau *is* dt:
    # this used to read complex(-dt.real,dt.imag), i.e. it negated the
    # real part only. Purely-imaginary dt (timeevolution.evolve_WF's
    # dt=1j*dt01, the only caller that existed before this was noticed)
    # was unaffected and is unchanged here; a *real* dt got exp(-dt*h),
    # so exponential()'s Hermitian branch computed e^{-h} where
    # edchain.exponential (algebra.expm(h)) and the user guide both say
    # e^{+h}. That was invisible for as long as the Hermitian branch was
    # unreachable for multi-site operators (see exponential() above), and
    # in examples/time_evolution/exponential_EV, whose sum(Sx) has the
    # same expectation value under either sign in a Z-polarized state.
    tau = complex(dt)
    handle = self._session.exponential_apply(h.to_terms(),wfa.cpp_handle,
            tau,int(nt0))
    return mps.MPS(self,cpp_handle=handle).copy()

# `if self.mode is not None: mode = self.mode` used to stand where
# resolve_mode() is called in the next two functions. That idiom honours
# an explicit sc.mode="ED" but is blind to mode.py's *automatic*
# fallbacks (no compiled extension for the requested C++ version;
# itensor_version=3 on a chain with fewer than 3 sites), which return
# "ED" without ever writing self.mode. On such a chain get_gs() hands
# back an ED State while these functions still took their DMRG branch,
# failing several frames deep with "'State' object has no attribute
# 'cpp_handle'". resolve_mode() -- rather than get_mode() -- is what
# these want: it sees both fallbacks, without get_mode()'s extra
# conserved-sector guard, which is about *which ground state answers*
# and not about an inner product between two states the caller already
# holds.
# These two take the wavefunction's own type, like every other
# primitive in this file (see wavefunction_mode above), and NOT
# resolve_mode(). Routing them through resolve_mode() instead is an
# infinite recursion, not merely the wrong backend: EDchain.overlap is
# `return wf1.dot(wf2)`, and mps.MPS.dot is `return self.MBO.overlap(...)`,
# so handing an MPS to the ED branch comes straight back here and goes
# round again. Reached in practice by NH-DMRG on a 2-site chain, where
# mode.py's automatic itensor_version=3 ns<3 fallback makes resolve_mode
# answer "ED" while nhdmrg -- which hand-rolls its own two-site sweep and
# never calls dmrg(), so the fallback does not apply to it -- is still
# holding genuine MPS objects with live cpp_handles
# (tests/test_nhdmrg_generalized.py::test_nhdmrg_generalized_v3_short_chain_does_not_crash,
# RecursionError). That is the general shape: an automatic fallback
# describes which solver computes a *ground state*, and says nothing
# about two wavefunctions the caller is already holding.
def overlap(self,wf1,wf2,mode=None):
    """Compute the overlap <wf1|wf2>"""
    mode = wavefunction_mode(wf1,mode=mode) # solver, see above
    if mode=="DMRG": return overlap_dmrg(self,wf1,wf2)
    return self.get_ED_obj().overlap(wf1,wf2)


def overlap_aMb(self,wf1,A,wf2,mode=None):
    """Compute the overlap <wf1|M|wf2>"""
    mode = wavefunction_mode(wf1,mode=mode) # solver, see above
    if mode=="DMRG": return overlap_aMb_dmrg(self,wf1,A,wf2)
    return wf1.dot(A*wf2) # workaround


def overlap_dmrg(self,wf1,wf2):
    """Compute the overlap between wavefunctions"""
    return self._session.overlap(wf1.cpp_handle,wf2.cpp_handle)


def overlap_aMb_dmrg(self,wf1,A,wf2):
    """Compute the overlap between wavefunctions"""
    from .multioperator import MultiOperator
    from .multioperatortk.staticoperator import StaticOperator
    if type(A)==StaticOperator:
        return A.aMb(wf1,wf2)
#        return wf1.dot(A*wf2) # workaround
    else:
        return overlap_aMb_dmrg_MO(self,wf1,A,wf2)



def overlap_aMb_dmrg_MO(self,wf1,A,wf2):
    """Compute the overlap between wavefunctions, with A a multioperator"""
    from .multioperator import obj2MO
    A = obj2MO(A) # convert to a MO
    return self._session.overlap_aMb(wf1.cpp_handle,A.to_terms(),wf2.cpp_handle)


def applyoperator(self,A,wf,mode=None,**kwargs):
    mode = wavefunction_mode(wf,mode=mode)
    if mode=="DMRG": return applyoperator_dmrg(self,A,wf)
    elif mode=="ED":
        return self.get_ED_obj().applyoperator(A,wf)


def applyinverse(self,A,wf,mode=None,**kwargs):
    mode = wavefunction_mode(wf,mode=mode)
    # note mode= is consumed above and deliberately *not* forwarded:
    # applyinverse_dmrg takes only delta/maxn, so passing it on was a
    # TypeError four frames deep
    if mode=="DMRG": return applyinverse_dmrg(self,A,wf,**kwargs)
    elif mode=="ED":
        return wf.applyinverse(A)
#        return self.get_ED_obj().applyoperator(A,wf)


def summps(self,wf1,wf2,mode=None,**kwargs):
    mode = wavefunction_mode(wf1,mode=mode)
    if mode=="DMRG": return summps_dmrg(self,wf1,wf2)
    elif mode=="ED": return wf1 + wf2 #self.get_ED_obj().summps(A,wf1,wf2)



def scale_mps(self,wf,x,mode=None):
    """Multiply an MPS by a number.

    Every backend whose session exposes scale_mps() (itensor_version 2, 3
    and "python") does this by rescaling a single site tensor -- O(chi^2 d),
    no contraction, no bond growth and no truncation. The fallback is the
    original formulation, "apply the operator x*Id", which is a full
    truncating MPO sweep over the whole chain just to multiply a
    wavefunction by a number; it is kept only for backends without the
    binding (mpsjulialive's session has no scale_mps, and cvm.py does run
    on julia_live).

    Why this matters enough to have its own primitive: measured with
    cProfile on a 16-site Heisenberg CVM run (cvm.py at cvm_maxm=40,
    itensor_version=3), the identity-MPO route accounted for 15.4 s of the
    43.8 s total -- 939 calls, i.e. ~5 per conjugate-gradient iteration
    (eta^2*v, alpha*p, alpha*Ap, the (-1)*x hidden inside MPS.__sub__, and
    beta*p). On itensor_version="python" the same measurement was 2.22 s of
    4.09 s. Note this is not only faster but strictly *less lossy*: the old
    route ran the MPS through a cutoff/maxdim compression on every scalar
    multiplication, so results shift slightly (within DMRG tolerance) when
    switching to this path.
    """
    if wavefunction_mode(wf,mode=mode)=="ED":
        # an ED State is a dense vector; scaling it needs no session at
        # all. Without this branch the State fell straight through to
        # wf.cpp_handle below and died with an AttributeError.
        return x*wf
    session = getattr(self,"_session",None)
    if session is not None and hasattr(session,"scale_mps"):
        handle = session.scale_mps(wf.cpp_handle,complex(x))
        return mps.MPS(self,cpp_handle=handle).copy()
    from .multioperator import identity
    return applyoperator(self,x*identity(),wf)


def summps_dmrg(self,wf1,wf2):
    """Apply operator to a many body wavefunction"""
    handle = self._session.sum_mps(wf1.cpp_handle,wf2.cpp_handle)
    return mps.MPS(self,cpp_handle=handle).copy()


def applyoperator_dmrg(self,A,wf):
    """Apply operator via the in-process pybind11 extension
    (mpscpp2/chain_session.h's Chain::apply_operator)."""
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_verbose(self.verbose)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    handle = self._session.apply_operator(A.to_terms(),wf.cpp_handle)
    return mps.MPS(self,cpp_handle=handle).copy()


def applyinverse_dmrg(self,A,wf,delta=None,maxn=None):
    """Apply operator to a many body wavefunction"""
    if delta is None: delta = self.cvm_tol # overwrite
    if maxn is None: maxn = self.cvm_nit # overwrite
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_verbose(self.verbose)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    handle = self._session.apply_inverse(A.to_terms(),wf.cpp_handle,
            delta,int(maxn))
    return mps.MPS(self,cpp_handle=handle).copy()



def operator_norm(self,op,ntries=5,simplify=True,mode=None):
    """Given a certain operator, compute its norm.

    mode= picks the solver the random witness states are drawn from, so
    this (and is_zero_operator on top of it) really does take the same
    mode= as the rest of the API rather than raising TypeError on it.
    The ED route works unchanged: State supports both op*wf and
    wf.overlap(wf)."""
    if simplify: op = op.simplify() # simplify the operator
    out = [] # empty list
    for i in range(ntries):
        if mode is None: wf = self.random_mps() # random wavefunction
        else: wf = self.random_mps(mode=mode)
        wf = op*wf # apply the operator
        o = (wf.overlap(wf)).real
        out.append(o)
    return np.mean(out) # return the norm



def is_hermitian(self,op):
    """Given a certain operator, check if it is Hermitian.

    This only needs to tell "dh := op-op.get_dagger() is exactly zero"
    apart from "dh is a genuine nonzero operator" -- it does not need its
    witness wavefunction to be numerically accurate, so building that
    witness (and applying dh to it) at the caller's *production* bond
    dimension self.maxm, as this used to do unconditionally, is pure
    waste: confirmed directly via cProfile, this check alone accounted
    for 31%-38% of gs_energy()'s own total wall time on representative
    chains (n=16/maxm=60: 0.275s of 0.881s; n=24/maxm=100: 1.108s of
    2.906s), dominated by random_mps()'s own MPS-sum compression sweep
    and applyoperator's MPO-application compression sweep, both bounded
    by self.maxm. A nonzero dh applied to *any* generic random state --
    even a low-bond-dimension one, down to a plain product state -- has
    overwhelming probability of nonzero norm: dh acts through the same
    fixed set of local terms regardless of the witness's bond dimension,
    and an exact cancellation on one specific low-entanglement witness
    state is a measure-zero coincidence for generic coefficients, no
    more or less likely than at self.maxm's own scale (this was already
    a single-random-sample probabilistic witness, not an exhaustive
    proof, even before this change). So the witness here is built and
    probed at a small, fixed bond dimension instead of self.maxm,
    restored immediately after via try/finally."""
    op = op - op.get_dagger()
    old_maxm = self.maxm
    self.maxm = min(old_maxm, 8)
    # ...and a bond dimension of 1 outright on a long chain: the witness is
    # built as a sum of two random states truncated back to self.maxm, and
    # every bond of that sum discards weight, so its surviving norm decays
    # roughly geometrically with the number of sites until it falls under
    # MPS.normalize()'s 1e-8 floor (which returns None -- see below). A
    # bond-dimension-1 witness truncates nothing and is explicitly enough
    # here, per this function's own docstring.
    if self.ns > 24: self.maxm = 1
    try:
        wf = self.random_mps() # random wavefunction, small bond dimension
        if wf is None:
            # random_mps() builds its witness as a *sum* of two random
            # states truncated to self.maxm, and MPS.normalize() returns
            # None (with a warning) below its own 1e-8 norm floor. Each
            # bond of that sum discards weight, so the surviving norm
            # falls off roughly geometrically with chain length, and past
            # ~40 sites it can land under the floor -- worst for a chain
            # in conserved-sector mode, whose start state is deliberately
            # a sum of near-orthogonal product states and so has no small
            # low-rank approximation at all. Confirmed directly: a
            # 40-site Heisenberg chain with set_conserved_sector(Sz=0)
            # crashed here every time, the same chain without a sector
            # only occasionally. A bond-dimension-1 witness is immune (no
            # truncation happens at all) and is explicitly enough for
            # this probe -- see this function's own docstring.
            self.maxm = 1
            wf = self.random_mps()
            if wf is None: return True # no usable witness: treat as Hermitian
        wf = op*wf # apply the operator
        # applyoperator() normalizes its result too, and a witness that
        # op-op^dagger annihilates is exactly what a Hermitian op looks
        # like here.
        if wf is None: return True
        norm = (wf.dot(wf)).real
    finally:
        self.maxm = old_maxm
    return not norm>1e-4






from .algebra.arnolditk import mpsarnoldi
from .algebra.arnolditk import lowest_energy as lowest_energy_arnoldi
from .algebra.arnolditk import lowest_energy_non_hermitian as lowest_energy_non_hermitian_arnoldi
from .algebra.arnolditk import gram_smith_single

from .algebra.arpacktk import mpsiram
from .algebra.arpacktk import lowest_energy as lowest_energy_iram
from .algebra.arpacktk import lowest_energy_non_hermitian as lowest_energy_non_hermitian_iram
from .algebra.arpacktk import excited_states as mps_excited_states
from .algebra.arpacktk import mpsiram_shift_invert
from .algebra.arpacktk import shift_invert_excited_states
from .algebra.arpacktk import mpsiram_generalized
from .algebra.arpacktk import generalized_excited_states

# IRAM (algebra/arpacktk.py, ported from ARPACK) is the default MPS
# Arnoldi solver: it reuses its compressed Krylov subspace across
# restarts instead of rebuilding it from scratch, needing fewer H|psi>
# applications than arnolditk's explicit-restart Arnoldi on most spectra
# -- see examples/non_hermitian/arnoldi_vs_iram_benchmark. The _arnoldi
# and _iram suffixed names above stay available for explicit selection or
# head-to-head comparison; these bare names are the new default.
lowest_energy = lowest_energy_iram
lowest_energy_non_hermitian = lowest_energy_non_hermitian_iram


def toMPO(self,H,mode="DMRG"):
    """Transport an operator into a matrix-product operator"""
    if mode=="DMRG":
        if self.itensor_version in (2,3,"python"):
            from .multioperatortk.staticoperator import StaticOperator
            return StaticOperator(H,self) 
        elif self.itensor_version=="julia_live":
            from .mpsjulialive.mpo import MPO
            return MPO(H,MBO=self)
        else: raise NotImplementedError("toMPO is not implemented for "
                "itensor_version="+repr(self.itensor_version))
    elif mode=="ED":
        from .edtk.edchain import EDOperator
        return EDOperator(H,self.get_ED_obj())
    else: raise ValueError("Unrecognized mode "+repr(mode))



def conjugate_mps(self,wf):
    """Apply operator to a many body wavefunction"""
    handle = self._session.conjugate(wf.cpp_handle)
    return mps.MPS(self,cpp_handle=handle).copy()


from .mpsalgebratk.trace import trace
from .mpsalgebratk.trace import inverse_trace


from .mpsalgebratk.disentangle import disentangle_manifold




