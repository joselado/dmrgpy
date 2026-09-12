

import numpy as np





def dynamical_correlator_cvm_explicit(self,name=None,
        delta=1e-1,es=np.linspace(0.,5.0,300)):
    """
    Compute the dynamical correlator by explicitly inverting (z-H) at
    z = omega+E0 +- i*delta, i.e. the most literal possible transcription
    of the correction-vector method -- which is what makes it the
    cross-check for a suspicious submode="CVM" curve.

    Returns dmrgpy's house convention (see src/dmrgpy/dynamics.py's
    module docstring), the complex Lehmann density

        i*(G^R - G^A)/(2*pi),   G^R/G^A = <GS|A (w+E0 +- i*delta-H)^-1 B|GS>

    so the same numbers submode="CVM"/"KPM"/"ROOTN"/"EX" return, as the
    user guide promises.
    """
    print("Non Hermitian mode in dynamical correlator")
    ### So far this just works for onsite correlators
    A,B = name[0],name[1]
    if not self.is_zero_operator(A.get_dagger()-B):
        # a bare `raise` here surfaced as "RuntimeError: No active
        # exception to reraise", naming neither the submode nor the
        # restriction it is complaining about
        raise NotImplementedError(
            "get_dynamical_correlator: submode='CVM_explicit' is only "
            "implemented for a Hermitian operator pair, A^dagger == B "
            "(it inverts (z-H) against the single vector B|GS>, so the "
            "two sides cannot differ). Use submode='CVM', 'KPM' or 'EX' "
            "for a general off-diagonal pair.")
    wf = self.get_gs() # get the ground state
    wfa = A.get_dagger()*wf # apply A to the GS
    wfb = B*wf # apply B to the GS
    e0 = self.gs_energy() # ground state energy
    Hp = self.hamiltonian - e0
    def f(e,delta): # function to compute
        wfi = self.applyinverse(-self.hamiltonian+(e0+e+1j*delta),wfa)
        return wfb.dot(wfi) # return result
#    from .analyticcontinuation import imag2real
    outz = np.array([f(z,delta) - f(z,-delta) for z in es]) # complex axis
    # outz = G^R - G^A, so the house convention is i*outz/(2*pi). This
    # used to return 1j*outz/np.pi, i.e. exactly TWICE the spectral
    # function it is documented to share with submode="CVM" -- the 1/2
    # that turns the two-sided resolvent difference into Im G was simply
    # missing, and no test ever compared the two submodes' amplitudes
    # (measured ratio 2.0001 against the exact ED reference on every
    # backend). The np.abs() that wrapped it is gone too: it is a no-op
    # on the Hermitian chains the guard above admits, and on a genuinely
    # non-Hermitian H it discarded the phase of a quantity every other
    # submode returns complex. (2026-09 audit, finding #11.)
    return es,0.5j*outz/np.pi




dynamical_correlator_non_hermitian = dynamical_correlator_cvm_explicit






