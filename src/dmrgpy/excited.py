import numpy as np
from . import mps
from scipy import linalg as lg


def get_excited_states_dmrg(self,n=2,noise=0.0,scale=10.0):
    """Return excited state energies"""
    self.get_gs()
    if self.excited_gram_schmidt: sm = "true"
    else: sm = "false"
    task = {"excited":"true",
            "nexcited":str(n),
            "noise":str(noise),
            "excited_gram_schmidt":sm,
            "scale_lagrange_excited":str(scale),
            }
    self.task = task
    self.write_task()
    self.write_hamiltonian() # write the Hamiltonian to a file
    self.run() # perform the calculation
    wfs = [] # read the wavefunctions
    for i in range(n):
        wf = mps.MPS(MBO=self,name="wavefunction_"+str(i)+".mps").copy() 
        wfs.append(wf) # store this one
    out = self.execute(lambda: np.genfromtxt("EXCITED.OUT").T)
    return out[0],wfs # return energies and wavefunctions 


def get_excited(*args,**kwargs):
    """Return the excited state energies"""
    (es,ws) = get_excited_states(*args,**kwargs)
    return es



def get_excited_states(self,n=2,purify=True,**kwargs):
    """Excited states"""
    H = self.hamiltonian # make a check
    if not self.is_hermitian(H): # non-Hermitian Hamiltonians
#        print("Non-Hermitian Hamiltonian, using Arnoldi method")
        return excited_states_non_hermitian(self,n=n,**kwargs)
    if n==1: # workaround for just the ground state
        w = self.get_gs(**kwargs)
        e0 = self.gs_energy(**kwargs)
        return ([e0],[w]) # return
    if not purify: # just compute excited states
        return get_excited_states_dmrg(self,n=n,**kwargs) # compute 
    else: # purify the states
        es,ws = get_excited_states_dmrg(self,n=n+2,**kwargs) # compute 
        ws = gram_smith(ws) # orthogonalize the MPS
        ne = len(es)
        h = np.zeros((ne,ne),dtype=np.complex)
        for i in range(ne):
          for j in range(ne):
              h[i,j] = ws[i].overlap(self.hamiltonian*ws[j])
        es = lg.eigvalsh(h) # redefine eigenvalues
        # redefine also the eigenvectors
        from .algebra.arnolditk import rediagonalize
        ws = rediagonalize(self.hamiltonian,ws) # rediagonalize
        return (es[0:n],ws[0:n])


from .algebra.arnolditk import gram_smith


def excited_states_non_hermitian(self,n=3,**kwargs):
    from . import mpsalgebra
    (es,wf) = mpsalgebra.mpsarnoldi(self,self.hamiltonian,mode="GS",
                nwf=n, # number of wavefunctions
                recursive_arnoldi=True,
                n = n + 2, # number of Krylov vectors
                **kwargs)
    if n==1: 
        es = [es]
        wf = [wf]  
    return (es,wf)


