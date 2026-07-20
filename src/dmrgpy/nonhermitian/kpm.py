"""Non-Hermitian Kernel Polynomial Method (NH-KPM) dynamical correlator.

Port of the biorthogonal Chebyshev-moment algorithm from NHKPM.jl
(https://github.com/GUANGZECHEN/NHKPM.jl, itself implementing the method of
Phys. Rev. Lett. 130, 100401) to dmrgpy's ED and mpscpp3 (ITensor v3)
backends. Unlike the ordinary (Hermitian) KPM dynamical correlator
(kpmdmrg.py), the ground state here is a biorthogonal right/left pair
(psi_R,psi_L) from the non-Hermitian solver (nhdmrg()/ED), and the
Chebyshev moments are recomputed from scratch at every requested frequency
(the expansion operator (z*Id-H)/E_max itself depends on z) rather than
amortized once over the whole spectrum as in the Hermitian case - this
mirrors the reference algorithm's own cost profile.
"""
import numpy as np
from ..algebra import kpm
from ..algebra import algebra


def dynamical_correlator_nhkpm(self,name=None,delta=1e-1,
        es=np.linspace(0.,5.0,300),E_max=None,n=200,kernel="jackson"):
    """NH-KPM dynamical correlator for the DMRG backends (mpscpp3 only for
    now). name=(A,B) computes <psi_L|A(z) B|psi_R> at z=e0+es+1j*delta,
    with (psi_L,psi_R) the biorthogonal ground state pair from nhdmrg().
    E_max (a real upper bound for the spectral radius of H) must be
    supplied by the user - there is no automatic estimator yet for
    non-Hermitian operators (see maximum_energy() in chain_session.h,
    which only works for Hermitian H)."""
    if name is None: raise ValueError("name=(A,B) must be provided")
    if E_max is None:
        raise ValueError("E_max (an upper bound for the spectral radius "
                "of the Hamiltonian) must be provided for the "
                "non-Hermitian KPM dynamical correlator")
    if self.itensor_version!=3:
        raise NotImplementedError("NH-KPM dynamical correlator is only "
                "implemented for itensor_version=3 so far, got "
                +str(self.itensor_version))
    from .. import multioperator
    e0 = self.gs_energy() # routes to nhdmrg, sets self.wf0/self.nh_left_wf
    psir = self.wf0
    psil = self.nh_left_wf
    # gs_energy_nhdmrg() renormalizes wf0=psir.normalize() to unit norm but
    # leaves nh_left_wf at nhdmrg()'s own <psil|psir>=1 scale, so the two
    # no longer satisfy <psil|psir>=1 together - restore that biorthogonal
    # convention here rather than assuming it still holds.
    norm = psil.dot(psir)
    psil = (1.0/np.conjugate(norm))*psil
    mi = name[1] # applied to the right ground state
    mj = name[0].get_dagger() # applied to the left ground state
    wfa = mi*psir
    wfb = mj*psil
    ident = multioperator.identity()
    out = np.zeros(len(es),dtype=np.complex128)
    for i in range(len(es)):
        z = e0+es[i]+1j*delta
        hs = (z*ident-self.hamiltonian)*(1./E_max)
        hs_dag = hs.get_dagger()
        mu_n = self._session.nhkpm_moments(hs.to_terms(),hs_dag.to_terms(),
                wfa.cpp_handle,wfb.cpp_handle,int(n),
                int(self.kpmmaxm),float(self.kpmcutoff))
        out[i] = kpm.spec_from_moments_nh(np.array(mu_n),kernel=kernel)
    return es,out


def dynamical_correlator_nhkpm_ed(self,name=None,delta=1e-1,
        es=np.linspace(0.,5.0,300),E_max=None,n=200,kernel="jackson"):
    """NH-KPM dynamical correlator for the ED backend. self is an EDchain
    (as passed in by edtk/dynamics.py). Computes <psi_L|A(z) B|psi_R> with
    (psi_L,psi_R) the biorthogonal ED ground state (algebra.
    biorthogonal_ground_state)."""
    if name is None: raise ValueError("name=(A,B) must be provided")
    if E_max is None:
        raise ValueError("E_max (an upper bound for the spectral radius "
                "of the Hamiltonian) must be provided for the "
                "non-Hermitian KPM dynamical correlator")
    from ..edtk.edchain import EDOperator
    A = EDOperator(name[0],self).SO
    B = EDOperator(name[1],self).SO
    h = self.get_hamiltonian()
    e0,vr,vl = algebra.biorthogonal_ground_state(h)
    wfa = B@vr
    wfb = algebra.dagger(A)@vl
    out = np.zeros(len(es),dtype=np.complex128)
    for i in range(len(es)):
        z = e0+es[i]+1j*delta
        out[i] = kpm.get_spec_kpm_nh(h,z,wfa,wfb,E_max,n=n,kernel=kernel)
    return es,out
