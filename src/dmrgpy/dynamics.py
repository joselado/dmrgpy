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
sum-rule (and delta->0) sense, which is what the test file pins it by,
and its `delta` is the same broadening: every KPM route picks its moment
count from algebra/kpm.py's polynomials_for_broadening, so the Jackson
line comes out at FWHM = 2*delta, the width the resolvent submodes give.
Two caveats, both measured rather than asserted. The width is exact at
the centre of the rescaled band and tightens as sqrt(1-x^2) towards its
edges, which is the kernel's own resolution profile and not a choice: no
single moment count gives one width across a whole band. For a
ground-state correlator that narrowing is the rule rather than the
exception, since E0 sits at x0 = -1/(2*kpm_scale) on every chain and the
band centre at E0 + W/2, so as W grows every intensive excitation tends to
0.70 of the requested width (0.157 at E0 on the kpm_energy_truncate
window); polynomials_for_broadening's docstring has the measured numbers
and the compensation, delta/sqrt(1-x(omega)^2). And the line is
Jackson-Gaussian, not Lorentzian, so at equal FWHM and equal integrated
weight its peak stands about 1.5x higher than the resolvent submodes'
(1.514, the kernel's own ratio on one pole at x=0); the 0.315 against
0.195 once quoted here, on a 6-site Heisenberg chain at delta=0.2, is
1.61 because its dominant pole sits at x=-0.527, where the line is
already 0.850 of 2*delta. What it is NOT any more is a
different curve on the two solvers: mode="ED" and mode="DMRG" used to
pick both the rescaling window and the moment count independently
(int(2*scale/delta) against round((emax-emin)/delta)*kpm_n_scale, on
windows differing by a factor of three), and disagreed pointwise by
4.4e-01 against a resolvent peak of 0.355 on a 4-site chain at
delta=0.15. They now share both on the default bandwidth-centred window
(mode="ED" refuses kpm_energy_truncate, which moves v3 and "python" onto
the ground-state-anchored one; 2026-09-24b audit, finding 4), and since
the DMRG routes were cut to
the calibrated n moments (they reconstructed from n+2; 2026-09-24 audit,
finding 4) they agree to 2.6e-04 there on itensor_version="python",
which is the linear interpolation of the moment reconstruction's own
10n-point grid onto es; on itensor_version=3 the median is the same and
the worst runs reach about 8e-3, set by the run-to-run noise of its
band-edge estimate emax. That was open item O2 of the 2026-09 audit. get_distribution()'s own KPM path
(kpmdmrg.general_kpm_moments) is deliberately NOT on this calibration:
it expands an arbitrary operator rather than the Hamiltonian, and its
delta still only sets a polynomial count.
submode="INV"/"CVM" under mode="ED" and submode="EX" have always computed
C_AB directly. The 2026-09 audit found four routes off the convention and
brought them onto it rather than the other way round: submode="ED",
submode="ROOTN" and mode="DMRG" submode="CVM" (finding #5), plus
submode="CVM_explicit", which returned exactly 2x C_AB on every backend
and additionally destroyed the sign of a negative-weight correlator with
an np.abs() (finding #11). submode="TD"/"TDZ", which returned the
complex one-sided Fourier transform, were brought onto it afterwards as
open item O1 of docs/audit_2026_09_hole_hunt.md; see the real-time
section below. The
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
(mpsjulialive/dynamics.py::_kpm_dynamical_correlator plus kpm.jl). Its
moment count is the shared one and, since the 2026-09-24 audit's finding
4, cut to exactly n like every other DMRG route; on the band-centre pole
of two decoupled dimers it matches mode="ED" to better than 1e-6 at the
reconstruction's own grid points (tests/test_audit_2026_09_24_kpm.py).
Nothing wider than that has been measured on it.

THE TWO REAL-TIME ROUTES: submode="TD" and submode="TDZ"
--------------------------------------------------------
These two are built on a one-sided transform, and a one-sided transform
is a resolvent rather than a density, so they used to be the exceptions
to everything above. They are not any more. A real-time run produces
C(t) for t>=0 only, and `_fourier_transform_correlator` turns that into
-(i/pi)*G^A_AB(w), whose real part is C_AB when Im M_n = 0 and whose
imaginary part is a dispersive term C_AB does not have, measured at 70%
of the density's own peak even on a Hermitian pair. The missing half of
the transform is not a second simulation: the backward half of the pair
(A,B) is the conjugate of the FORWARD half of the pair
(B^dagger,A^dagger), so

    C_AB = ( F[(A,B)] + conj(F[(B^dagger,A^dagger)]) ) / 2

with F the one-sided transform, and that collapses to Re F exactly when
A is provably B^dagger, which is every example in the documentation and
costs one evolution rather than two. See
timedependent.lehmann_density_from_one_sided, which both submodes share;
canonical.is_dagger_pair is the test, and it refuses rather than guesses
for a name with no known adjoint, which costs a second evolution; that
second evolution is built from get_dagger(), so it is right exactly when
get_dagger() knows the name's adjoint (it did not for ISy until the
2026-09-24 audit's finding 2 gave get_dagger() a phase).

Measured against an exact Lehmann sum built by dense diagonalization
outside dmrgpy, on the 2026-09 audit's own seeded 4-site complex-hopping
chain (A = Cdag_0, B = C_2, max|Im M_n| = 0.271, exact peak 0.2313,
delta=0.4, dt=0.1), max|y - exact|:

    submode="TD",  mode="DMRG"   1.16e-01  ->  2.57e-04  ->  3.2e-05
    submode="TD",  mode="ED"     2.76e-01  ->  2.57e-04  ->  3.1e-05
    submode="TDZ", mode="DMRG"   1.13e-01  ->  3.31e-02  ->  5.4e-04

the last column being after the 2026-09-24 audit's finding 7, below,

and on a 4-site Heisenberg chain with the Hermitian pair A = B = Sz_0
(peak 0.1869, delta=0.3) the imaginary part goes from 1.11e-01, 60% of
that peak, to exactly zero, with max|y - exact| going 1.11e-01 ->
2.82e-04. The 3.31e-02 TDZ kept after O1 was not its contour, as this
docstring used to say, but the frequency stage's linear interpolation of
an FFT grid of spacing 2*pi/(nt*dt) = 1.05*delta, shared with TD at
predict=False; since the damped sum is evaluated at each requested
frequency (2026-09-24 audit, finding 7) TDZ and TD at predict=False sit
at the 5.4e-04 finite-window floor, and the contour's own share,
max|y_TDZ - y_TD(predict=False)|, is about 1e-06. sxt_to_skomega, left on
the FFT stage then, came onto the direct sum with the 2026-09-24b audit's
finding 7. The mode="ED" row moved further than the others
because that route additionally read the operator pair in the opposite
order, so the two solvers computed the correlator of different pairs
under one submode name, invisible whenever A and B are the same
operator; see edtk/timedependent.evolution_DC.

This closes open item O1 of docs/audit_2026_09_hole_hunt.md. Numbers
change accordingly: any submode="TD"/"TDZ" spectrum from before this is
not comparable, its imaginary part most of all. The infinite-chain
S(k,omega) route (timedependent.sxt_to_skomega, used by
pyitensor.idmrg_window and infinitechain.py) still returns the raw
transform, because the reduction it performs never sees the operator
pair. It is on the house sign, though: its S(x,t) carries e^{-i(H-E_0)t}
and it conjugates the momentum series after the spatial sum, so its lines
sit at omega = +D_n and its real part is the house density whenever the
momentum-resolved weights are real; until the 2026-09-24b audit's
finding 8 every infinite-chain S(k,omega) was mirrored in omega.
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
        # the chain's current ground state, on the Python side and on the
        # session alike, and the Hamiltonian on the session: free on a
        # solved chain, the hand-off of a state set_gs()/set_initial_wf()
        # injected, a solve otherwise (groundstate.ground_state_on_session)
        from .groundstate import ground_state_on_session
        hermitian = self.is_hermitian(self.hamiltonian)
        # The KPM route's own argument checks go ahead of the ground state,
        # and SECTOR, which solves both its sectors on a clone and never
        # reads the caller's state, makes no solve here at all: with the
        # solve first, a malformed KPM call paid a full ground-state solve
        # before raising, and the first SECTOR call one it never read
        # (2026-09-24c audit, finding 11).
        if submode=="KPM" and hermitian: kpmdmrg.check_kpm_call(self,**kwargs)
        if submode!="SECTOR": ground_state_on_session(self)
        if not hermitian: # non Hermitian Hamiltonian
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



