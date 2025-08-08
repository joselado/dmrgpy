import numpy as np

# library to compute degeneracies of eigenstates


def gs_degeneracy(self,mode="DMRG",**kwargs):
    """Compute the degeneracy of the ground state"""
    if mode=="DMRG": # for DMRG, it depends
        if self.is_hermitian(self.hamiltonian): # brute force should be ok
            return gs_degeneracy_simple(self,mode=mode,**kwargs)
        else: # non Hermitian is more tricky (due to the excited states)
            e0 = self.gs_energy(**kwargs) # ground state energy
            return eigenvalue_degeneracy(self,self.hamiltonian,1.2*e0,**kwargs)
    elif mode=="ED": # brute force should be fine
        return gs_degeneracy_simple(self,mode=mode,**kwargs)



def gs_degeneracy_simple(self,dmode="real",delta=1e-2,n=1,**kwargs):
    """Compute the degeneracy of the ground state with a minimal algorithm"""
    while True:
        es = self.get_excited(n=n+2,**kwargs) # get the eigenenergies
        if dmode=="real":
            emin = np.min(es.real) # ground state energy
            des = np.abs(es.real-emin)**2 # distance to the minimum energy
        elif dmode=="complex": # for complex energies
            esre = es.real # real part of the energy
            emin = np.min(esre) # ground state energy (real part)
            ind = np.where(esre==emin)[0][0]# index in the array
            ezmin = es[ind] # complex energy of the ground state
            des = np.abs(es-ezmin)**2 # distance to the minimum energy
        else: raise
        deg = np.sum(np.exp(-(des/delta)**2)) # return degeneracy
        if deg<n: break # consider it done
        n = n + 1 # add one
    return deg # return the degeneracy





def eigenvalue_degeneracy(self,A,e,n=3,emode="real",delta=1e-2):
    """Given an operator and a certain eigenvalue, estimate
    what is the degeneracy"""
    from .algebra.arnolditk import mpsarnoldi
    while True: # this is mostly a DMRG implementation
        es,ws = mpsarnoldi(self,A,e=e,delta=delta,mode="ShiftInv",
            recursive_arnoldi=False,nwf=n+2,n0=0,
            verbose=1,maxit=20,nkry_min = 2*n
            )
        if emode=="real":
            emin = np.min(es.real) # ground state energy
            des = np.abs(es.real-emin)**2 # distance to the minimum energy
        else: raise
        deg = np.sum(np.exp(-(des/delta)**2)) # return degeneracy
        if deg<n: break # consider it done
        n = n + 1 # add one
    return deg







def pole_eigenvalue_degeneracy(self,A,e,delta=1e-1,**kwargs):
    """Compute the degeneracy of an eigenvalue using Cauchy's theorem"""
    # delta controls the "radious" of the area in the complex plane
    B = A - e # shift to the origin, so that the poles are there
    thetas = np.linspace(0.,1.,20) # angles around the origin
    zs = delta*np.exp(1j*2*np.pi*thetas) # complex number
    fs = [self.inverse_trace(B-z,delta=1e-3,**kwargs) for z in zs] # function
    ldim = self.get_logdimension() # return the dimensionality of the problem
    fs = np.array(fs) # convert to array
    intz = np.sum(fs*zs)/len(fs) # integral
    pole = -intz # compute pole
    nt = len(zs) # number of thetas
    dzs = [zs[(i+1)%nt] - zs[(i-1)%nt] for i in range(len(zs))] # finite differences
    pole2 = 1j*np.sum(fs*dzs)/4./np.pi # integral, but now with finite difference
#    print(pole,pole2)
    #pole = (pole + pole2)*np.exp(ldim)/2. # return the value of the pole
    pole = pole*np.exp(ldim) # return the value of the pole
    return pole





