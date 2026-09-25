

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
    user guide promises, for any pair (A,B): (z-H) is inverted against
    B|GS> and the result dotted with A^dagger|GS>, which is <GS|A.

    It used to invert against A^dagger|GS> and dot with B|GS>, i.e. to
    compute <GS|B^dagger (z-H)^-1 A^dagger|GS>, the density of the adjoint
    pair (B^dagger,A^dagger), and so it refused every pair but A^dagger ==
    B, where the two coincide. That gate was an absolute test on a squared
    norm (Many_Body_Chain.is_zero_operator's 1e-4, EDchain's 1e-8 on
    Tr(D D^dagger)), which admitted any pair of small enough operators and
    returned the adjoint pair's curve for it: 2.004 of the peak off for
    (1e-2*Sx0, 1e-2*Sy1) on every DMRG backend, 2.064 for (1e-5*Sx0,
    1e-5*Sy1) on a 6-site non-Hermitian ED chain (2026-09-25b audit,
    finding 16). With the two vectors the right way round the formula
    holds for every pair, so there is no gate left to get wrong; for an
    adjoint pair the two vectors are the same state and nothing changed.
    On a non-Hermitian H, <GS| is still the conjugate of the same right
    eigenvector on both sides, as it always was here.
    """
    print("Non Hermitian mode in dynamical correlator")
    A,B = name[0],name[1]
    wf = self.get_gs() # get the ground state
    wfa = A.get_dagger()*wf # A^dagger|GS>, whose conjugate is <GS|A
    wfb = B*wf # B|GS>
    e0 = self.gs_energy() # ground state energy
    def f(e,delta): # <GS|A (z-H)^-1 B|GS>, z = w+E0+i*delta
        wfi = self.applyinverse(-self.hamiltonian+(e0+e+1j*delta),wfb)
        return wfa.dot(wfi) # return result
#    from .analyticcontinuation import imag2real
    outz = np.array([f(z,delta) - f(z,-delta) for z in es]) # complex axis
    # outz = G^R - G^A, so the house convention is i*outz/(2*pi). This
    # used to return 1j*outz/np.pi, i.e. exactly TWICE the spectral
    # function it is documented to share with submode="CVM" -- the 1/2
    # that turns the two-sided resolvent difference into Im G was simply
    # missing, and no test ever compared the two submodes' amplitudes
    # (measured ratio 2.0001 against the exact ED reference on every
    # backend). The np.abs() that wrapped it is gone too: it was a no-op
    # on a Hermitian chain with an adjoint pair, the only case the old
    # gate admitted, and on a genuinely non-Hermitian H it discarded the
    # phase of a quantity every other submode returns complex. (2026-09
    # audit, finding #11.)
    return es,0.5j*outz/np.pi




dynamical_correlator_non_hermitian = dynamical_correlator_cvm_explicit






