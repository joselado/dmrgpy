"""Dispatch of get_dynamical_correlator over its submodes.

THE HOUSE CONVENTION
--------------------
The quantity this codebase calls a dynamical correlator -- what every
submode except "TD"/"TDZ" returns, on whichever of mode="DMRG"/mode="ED"
implements it -- is the *complex* Lehmann density of the operator pair
(A,B), Lorentzian-broadened by `delta`:

    C_AB(w) = sum_n M_n * delta/(pi*((w-D_n)^2+delta^2))
            -> sum_n M_n delta(w-D_n)   as delta->0
            =  i*(G^R_AB(w) - G^A_AB(w))/(2*pi)

with M_n = <GS|A|n><n|B|GS> and D_n = E_n-E_0, and

    G^R_AB(w) = <GS| A (w+E_0+i*delta-H)^{-1} B |GS>,  G^A = the same at
    -i*delta.

Note M_n is complex in general, so C_AB is complex too; it is real
exactly when every M_n is real. A Hermitian pair A = B^dagger guarantees
that (there M_n = |<n|B|GS>|^2, so C_AB is real *and* non-negative), and
it is the overwhelmingly common case, the one every example and most
tests use -- but it is not the only way: a real Hamiltonian with real
operators has real M_n for any pair, Hermitian or not. `Im M_n == 0`,
not `A == B^dagger`, is the discriminant wherever this docstring says two
conventions coincide.

The reason this, and not the also-common -(1/pi) Im G^R_AB, is the
convention here: -(1/pi) Im G^R = sum_n [Re(M_n)*delta -
Im(M_n)*(w-D_n)]/(pi*D) coincides with C_AB only for real M_n, and its
dispersive second term has principal-value tails that leak arbitrarily
far outside any finite frequency window, so it does not satisfy the
kernel-independent sum rule

    integral dw C_AB(w) = sum_n M_n = <GS|A B|GS>

that every submode on this convention can be checked against (see
tests/test_audit_2026_09_correlator-conventions.py). The default
submode="KPM" -- a Chebyshev expansion of the spectral density, which has
no notion of a retarded resolvent at all -- is on this convention in the
sum-rule (and delta->0) sense, which is what the test file pins it by; it
is NOT pointwise interchangeable with the resolvent submodes at a given
delta, since delta there only sets the polynomial count (edtk/dynamics.py
picks npol from int(2*scale/delta), and the DMRG side picks its own), so
measured peak heights differ by a factor of 2-3 at delta=0.15..0.6.
submode="INV"/"CVM" under mode="ED" and submode="EX" have always computed
C_AB directly. The 2026-09 audit found four routes off the convention and
brought them onto it rather than the other way round: submode="ED",
submode="ROOTN" and mode="DMRG" submode="CVM" (finding #5), plus
submode="CVM_explicit", which returned exactly 2x C_AB on every backend
and additionally destroyed the sign of a negative-weight correlator with
an np.abs() (finding #11). submode="TD"/"TDZ" are the routes still off
this convention -- they return the complex one-sided Fourier transform,
whose REAL part is C_AB when Im M_n == 0 -- recorded as open item O1 in
docs/audit_2026_09_hole_hunt.md rather than changed in that pass. The
test file above pins ED/{ED,INV,CVM,ROOTN} and DMRG/{CVM,ROOTN}
pointwise against an exact Lehmann sum, and ED/{ED,INV,CVM,ROOTN,KPM}
plus DMRG/KPM against the sum rule.

The two backends do not offer the same submodes, so "every submode"
never means "every (mode,submode) pair": mode="DMRG" takes the names in
SUBMODES below, while mode="ED" implements KPM, ED, EX, INV, CVM, ROOTN
and TD and raises NotImplementedError for the rest (edtk/dynamics.py).
"INV" is ED-only; "TDZ", "CVM_explicit", "CVMimag", "SECTOR" and
"maxent" are DMRG-only.

itensor_version="julia_live" is a third case and is only partly covered
by the statement above. Its CVM/TDZ/EX/maxent go through the very same
shared modules as every other backend (cvm.py, tdz.py, dcex.py,
distribution.py -- so they inherit whatever those return, the 2026-09
fixes included), but its KPM is its own implementation
(mpsjulialive/dynamics.py::_kpm_dynamical_correlator plus kpm.jl), which
that audit did not exercise: it mirrors kpmdmrg.py and shares the same
moment reconstruction, but nothing here has measured it, so this
docstring makes no claim about it.

THE TWO EXCEPTIONS: submode="TD" and submode="TDZ"
-------------------------------------------------
These two do NOT return C_AB, and finding #5's convention decision
explicitly did not touch them. Both end in
`timedependent._fourier_transform_correlator`, which returns the full
*complex* one-sided Fourier transform of the real-time correlator --
in the long-time limit -(i/pi)*G^A_AB(w), whose real part is C_AB when
Im M_n = 0 and whose imaginary part is the dispersive -(1/pi)*Re G^A_AB
that C_AB does not have. Measured on a 6-site Heisenberg chain with the
Hermitian pair A = B = Sz_0 (so M_n is real and C_AB is the ordinary
real density, peak 0.1421, delta=0.3): Re y reproduces C_AB to 2e-4 and
y reproduces -(i/pi)*G^A to 2e-4, while max|Im y| = 0.0996, i.e. 70% of
the density's own peak (submode="TDZ": 68%, with its complex-time
contour putting Re y 1.8e-2 from C_AB). So take `np.real(...)` of a
TD/TDZ result before comparing it against any other submode, and do not
read its imaginary part as a complex Lehmann weight. This is recorded as
an open item in docs/audit_2026_09_hole_hunt.md: moving them onto the
convention changes numbers on the most commonly used real-time route,
and was deliberately not done in that pass.
"""
from . import kpmdmrg
from . import timedependent
from . import cvm
from . import dcex
from . import tdz
from . import rootndmrg

SUBMODES = ("KPM","TD","TDZ","CVM","CVM_explicit","CVMimag","ROOTN","EX",
            "SECTOR","maxent")


def get_dynamical_correlator(self,submode="KPM",**kwargs):
    if submode not in SUBMODES:
        # checked up front: otherwise an unrecognized submode on a
        # non-Hermitian Hamiltonian gets described as "Hermitian-only"
        # below, and on a Hermitian one it used to reach a bare `raise`
        raise ValueError(
            "get_dynamical_correlator: unrecognized submode "+repr(submode)
            +"; expected one of "+", ".join(SUBMODES))
    if self.itensor_version in (2,3,"python"): # C++ or pure-Python
        self.set_initial_wf(self.wf0) # set the initial wavefunction
        if not self.is_hermitian(self.hamiltonian): # non Hermitian Hamiltonian
            # Per-submode, not wholesale. This check used to run before the
            # dispatch below and return the explicit resolvent for
            # *everything* except "KPM", so on a non-Hermitian Hamiltonian
            # the caller's submode= was a no-op: "EX", "maxent", "ROOTN",
            # "TD", "CVM" and "CVM_explicit" all returned the same curve
            # (bit-identical on the ED path). The same substitution reached
            # mode="ED", submode="ED" -- i.e. the exact Lehmann reference a
            # user would cross-validate against quietly became the
            # approximate resolvent and agreed with itself. The julia_live
            # branch below already fixed this shape for its own backend;
            # this is the same fix for (2,3,"python"), see also
            # edtk/dynamics.py.
            if submode=="KPM": # non-Hermitian KPM
                from .nonhermitian.kpm import dynamical_correlator_nhkpm
                return dynamical_correlator_nhkpm(self,**kwargs)
            elif submode in ("CVM","CVM_explicit"):
                # not a substitution: this *is* the non-Hermitian
                # implementation of the correction-vector resolvent both
                # of these compute
                from .nonhermitian.dynamics import dynamical_correlator_non_hermitian
                return dynamical_correlator_non_hermitian(self,**kwargs)
            elif submode in ("EX","maxent"):
                pass # backend-agnostic, non-Hermitian-capable: fall through
            else:
                raise NotImplementedError(
                    "get_dynamical_correlator: submode=%r assumes a "
                    "Hermitian Hamiltonian and has no non-Hermitian "
                    "implementation (KPM/CVM/CVM_explicit do, and EX/maxent "
                    "are non-Hermitian-capable already). It used to return "
                    "the CVM_explicit resolvent instead, silently."%submode)
        if submode=="KPM": # KPM method
            return kpmdmrg.get_dynamical_correlator(self,**kwargs)
        elif submode=="TD": # time dependent
            return timedependent.dynamical_correlator(self,**kwargs)
        elif submode=="TDZ": # complex-time evolution (arXiv:2311.10909)
            return tdz.dynamical_correlator_tdz(self,**kwargs)
        elif submode=="CVM_explicit": # CVM mode
            return cvm.dynamical_correlator_cvm_explicit(self,**kwargs)
        elif submode=="CVM": # CVM mode
            return cvm.dynamical_correlator(self,**kwargs)
        elif submode=="CVMimag": # CVM mode
            return cvm.dynamical_correlator_analytic_continuation(self,**kwargs)
        elif submode=="ROOTN": # root-N Krylov correction-vector
            return rootndmrg.dynamical_correlator(self,**kwargs)
        elif submode=="EX": # EX mode
            return dcex.dynamical_correlator(self,**kwargs)
        elif submode=="SECTOR": # sector-resolved Lehmann sum
            # itensor_version=2 reaches here (this branch admits 2, 3 and
            # "python") but has no quantum numbers at all; sectordc's own
            # _check_backend says so, naming the two backends that do.
            from . import sectordc
            return sectordc.dynamical_correlator(self,**kwargs)
        elif submode=="maxent": # Max ent mode
            from .distribution import dynamical_correlator_positive_defined
            return dynamical_correlator_positive_defined(self,**kwargs)
        else:
            # a bare `raise` here used to surface as the thoroughly
            # unhelpful "RuntimeError: No active exception to reraise"
            raise ValueError(
                "get_dynamical_correlator: unrecognized submode %r; "
                "expected one of KPM, TD, TDZ, CVM, CVM_explicit, CVMimag, "
                "ROOTN, EX, SECTOR, maxent"%submode)
    elif self.itensor_version=="julia_live": # Julia version
        # Only KPM/CVM/TDZ assume a Hermitian Hamiltonian (Chebyshev
        # spectrum-in-[-1,1], resolvent CG solve, and the TDZ damping
        # mechanism respectively) -- EX and maxent are already
        # backend-agnostic MultiOperator/MPS algebra with their own
        # working non-Hermitian path (dcex.py -> excited.py's
        # excited_states_non_hermitian, not itensor_version-gated) and
        # must not be blocked here. This check used to run before submode
        # dispatch entirely, which also rejected EX/maxent even though
        # they work fine -- confirmed directly, this was blocking a
        # working code path.
        if submode=="SECTOR":
            # No quantum numbers on this backend at all, so no sectors:
            # say that, rather than letting mpsjulialive fail on an
            # unrecognized submode several frames deeper.
            raise NotImplementedError(
                "get_dynamical_correlator: submode=\"SECTOR\" needs "
                "conserved-sector support, which itensor_version="
                "'julia_live' does not have -- use itensor_version=3 or "
                "itensor_version=\"python\"")
        if submode in ("KPM","CVM","TDZ") and not self.is_hermitian(self.hamiltonian):
            # unlike the (2,3,"python") branch above, there is no
            # non-Hermitian route to fall back to here for these
            # submodes: dynamical_correlator_non_hermitian ultimately
            # needs applyinverse_dmrg(), which is self._session-only
            # (mpsalgebra.py) and also dispatches on type(wf)==mps.MPS --
            # the *top-level* MPS class, not mpsjulialive.mps.MPS -- so it
            # would fail regardless. Silently running the Hermitian-only
            # KPM/CVM/TDZ math on a non-Hermitian Hamiltonian produces
            # numerically wrong output with no error; raise instead.
            raise NotImplementedError(
                "get_dynamical_correlator: itensor_version='julia_live' "
                "does not implement non-Hermitian Hamiltonians for "
                "submode=%r (KPM/CVM/TDZ all assume a Hermitian one); use "
                "submode='EX'/'maxent', or itensor_version in "
                "(2,3,'python') instead"%submode)
        from .mpsjulialive import dynamics as dynamicsjl
        return dynamicsjl.get_dynamical_correlator(self,submode=submode,**kwargs)
    else: raise



