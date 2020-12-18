from .algebra.kpm import generate_profile
import numpy as np



def get_moments_distribution(self,X=None,num_p=None):
  """Get the mments with DMRG"""
  # define the dictionary
  self.gs_energy()
  task = {      "distribution": "true",
                "kpmmaxm":str(self.kpmmaxm),
                "kpm_scale":str(self.kpm_scale),
                "kpm_accelerate":self.kpm_accelerate,
                "kpm_n_scale":str(self.kpm_n_scale),
                "kpm_num_polynomials":str(num_p),
                "kpm_cutoff":str(self.kpmcutoff),
                }
  self.execute(lambda: X.write(name="kpm_distribution_multioperator.in")) 
  self.task = task # assign tasks
  self.write_task()
  self.write_hamiltonian() # write the Hamiltonian to a file
  self.run() # perform the calculation
  m = self.execute(lambda: np.genfromtxt("KPM_MOMENTS.OUT").transpose())
  return m[0]+1j*m[1]





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
    try: from .maxenttk.pymaxent import reconstruct
    except:
        print("Not functional yet")
        exit()
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
    if not (A-B.get_dagger()).is_zero(): 
        print("Non positive correlator")
        raise
    e0 = self.gs_energy() # ground state energy
    h0 = self.hamiltonian-e0*(1.+0j) # shift Hamiltonian
    wf = self.get_gs() # get ground state
    return get_distribution_maxent(self,X=h0,wf=B*wf,
            bnds=[-0.1,5.0],x=es,n=6)

