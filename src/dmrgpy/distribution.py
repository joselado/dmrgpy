from .algebra.kpm import generate_profile
import numpy as np



def get_distribution(self,**kwargs):
    """
    Compute a dynamical correlator using the KPM-DMRG method
    """
    from .kpmdmrg import general_kpm
    return general_kpm(self,**kwargs)



def get_distribution_moments(self,**kwargs):
    """
    Compute a dynamical correlator using the KPM-DMRG method
    """
    from .kpmdmrg import general_kpm_moments
    return general_kpm_moments(self,**kwargs)



def get_distribution_maxent(self,X=None,wf=None,
        n=10,x=None,bnds=[-3.0,3.],**kwargs):
    """
    Compute a distribuion using a maxentropy method
    """
    if wf is None: wf = self.get_gs(**kwargs) # get wavefunction
    from .maxenttk.pymaxent import reconstruct
    from .vev import power_vev
    mu = power_vev(self,n=n,X=X,wf=wf).real
    scale = mu[0] # scale of the problem (1 for probabilities)
    mu = mu/scale # normalize
#    mu = [self.vev(X,npow=i).real for i in range(n)] # compute moments
    sol, lambdas = reconstruct(mu,bnds=bnds)
    if x is None: x = np.linspace(-1.,1.,300)
    return x,sol(x)*scale


def dynamical_correlator_positive_defined(self,name=None,
        es=np.linspace(-1.0,5.0,400),**kwargs):
    """Return a dynamical correlator that is positive defined"""
    A,B = name[0],name[1]
    # canonical.is_dagger_pair rather than (A-B.get_dagger()).is_zero():
    # get_dagger() leaves an operator name it does not recognize exactly
    # as it found it, so such an operator cancels against its own
    # "dagger" and passes this test whatever it is, which is the one
    # thing this test exists to catch. is_dagger_pair refuses instead
    # when either side names something whose adjoint is unknown.
    from .multioperatortk import canonical
    if not canonical.is_dagger_pair(A,B):
        raise ValueError("dynamical_correlator_positive_defined needs "
                "name=(A,B) with A = B^dagger, which could not be "
                "established for this pair: the distribution it "
                "reconstructs is only positive definite then")
    e0 = self.gs_energy() # ground state energy
    h0 = self.hamiltonian-e0*(1.+0j) # shift Hamiltonian
    wf = self.get_gs() # get ground state
    return get_distribution_maxent(self,X=h0,wf=B*wf,
            bnds=[-0.1,5.0],x=es,n=6)

