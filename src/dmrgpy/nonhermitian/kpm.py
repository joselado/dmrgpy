"""Non-Hermitian Kernel Polynomial Method (NH-KPM) dynamical correlator.

Port of the biorthogonal Chebyshev-moment algorithm from NHKPM.jl
(https://github.com/GUANGZECHEN/NHKPM.jl, itself implementing the method of
Phys. Rev. Lett. 130, 100401) to dmrgpy's ED and mpscpp3 (ITensor v3)
backends, and to "python" and julia_live. Unlike the ordinary
(Hermitian) KPM dynamical correlator (kpmdmrg.py), the ground state here
is a biorthogonal right/left pair (psi_R,psi_L) from the non-Hermitian
solver (nhdmrg()/ED), and the Chebyshev moments are recomputed from scratch at every requested frequency
(the expansion operator (z*Id-H)/E_max itself depends on z) rather than
amortized once over the whole spectrum as in the Hermitian case - this
mirrors the reference algorithm's own cost profile.
"""
import numpy as np
from ..algebra import kpm
from ..algebra import algebra


def dynamical_correlator_nhkpm(self,name=None,delta=1e-1,
        es=np.linspace(0.,5.0,300),E_max=None,n=200,kernel="jackson"):
    """NH-KPM dynamical correlator for the DMRG backends (itensor_version
    3, "python" and "julia_live"). name=(A,B) computes <psi_L|A(z) B|psi_R> at
    z=e0+es+1j*delta, with (psi_L,psi_R) the biorthogonal ground state
    pair from nhdmrg(). E_max (a real upper bound for the spectral radius
    of H) must be supplied by the user - there is no automatic estimator
    yet for non-Hermitian operators (see maximum_energy() in
    chain_session.h, which only works for Hermitian H). A state set by
    hand has no left partner, so after set_gs()/set_initial_wf()/
    gs_energy(wf0=...) this raises rather than pairing it with the left
    state of an earlier solve."""
    if name is None: raise ValueError("name=(A,B) must be provided")
    if E_max is None:
        raise ValueError("E_max (an upper bound for the spectral radius "
                "of the Hamiltonian) must be provided for the "
                "non-Hermitian KPM dynamical correlator")
    if self.itensor_version not in (3,"python","julia_live"):
        raise NotImplementedError("NH-KPM dynamical correlator is only "
                "implemented for itensor_version 3, \"python\" or "
                "\"julia_live\" so far, got "+str(self.itensor_version))
    from .. import multioperator
    # routes to NH-DMRG and sets self.wf0/self.nh_left_wf, or returns the
    # stored pair; on julia_live through mpsjulialive/groundstate.py
    e0 = self.gs_energy()
    psir = self.wf0
    psil = getattr(self,"nh_left_wf",None)
    if psil is None or getattr(self,"_nh_left_for",None) is not psir:
        # A right state with no left state of its own: set_gs(),
        # set_initial_wf() and gs_energy(wf0=x, reconverge=False) give only
        # the right one, and a right state that is not an eigenstate has no
        # biorthogonal partner at all. This used to pair x with the
        # nh_left_wf of the last NH-DMRG solve, the left eigenvector of a
        # different state, and renormalize by <psil|x> (0.93 on a 4-site
        # chain), or raise AttributeError when there had been no solve.
        raise RuntimeError("the non-Hermitian KPM correlator needs the "
                "left eigenvector that pairs with the chain's right state, "
                "and this chain has none: its state was set with set_gs(), "
                "set_initial_wf() or gs_energy(wf0=...), which give a right "
                "state only. Solve the biorthogonal pair with NH-DMRG "
                "instead, restart() and then gs_energy(). On julia_live a "
                "solve whose left and right vectors came for two "
                "eigenvalues tied in the real part keeps no left one "
                "either; break the tie in H")
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
        if self.itensor_version=="julia_live":
            # mpsjulialive/nhkpm.jl, the same recursion on the Julia MPOs
            # of hs and of its dagger, built at the MultiOperator level
            from ..mpsjulialive.mpo import get_MPO
            from ..mpsjulialive.juliasession import Main
            mu_n = Main.nhkpm_moments(get_MPO(hs,self),get_MPO(hs_dag,self),
                    wfa.jlmps,wfb.jlmps,int(n),
                    int(self.kpmmaxm),float(self.kpmcutoff))
            mu_n = [complex(x) for x in mu_n]
        else:
            mu_n = self._session.nhkpm_moments(hs.to_terms(),
                    hs_dag.to_terms(),wfa.cpp_handle,wfb.cpp_handle,int(n),
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
