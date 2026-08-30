from . import kpmdmrg
from . import timedependent
from . import cvm
from . import dcex
from . import tdz
from . import rootndmrg

SUBMODES = ("KPM","TD","TDZ","CVM","CVM_explicit","CVMimag","ROOTN","EX",
            "maxent")


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
        elif submode=="maxent": # Max ent mode
            from .distribution import dynamical_correlator_positive_defined
            return dynamical_correlator_positive_defined(self,**kwargs)
        else:
            # a bare `raise` here used to surface as the thoroughly
            # unhelpful "RuntimeError: No active exception to reraise"
            raise ValueError(
                "get_dynamical_correlator: unrecognized submode %r; "
                "expected one of KPM, TD, TDZ, CVM, CVM_explicit, CVMimag, "
                "ROOTN, EX, maxent"%submode)
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



