# routines to perform mean field calculations
import numpy as np


def spinchain_meanfield(sc,p=0.0,mix=0.9,m0=None,maxerror=1e-06,
        mode="default"):
    """Mean field calculation of a spin chain,
    p controls the degree of many body correlations"""
    sc0 = sc.copy() # copy the spin chain object
    if m0 is None:
        fold = np.array([np.random.random(3) for i in range(sc.ns)]) 
    else: fold = np.array(m0)
    def get_new(mold):
        """Return new magnetization"""
        sc = sc0.copy() # copy the initial Hamiltonian
        def m2mf(mnew):
          fnew = np.array([np.zeros(3) for i in range(sc.ns)]) # new mean field
          for c in sc.exchange: # loop over coupled sites
              fnew[c.i] += Av(c.g,mnew[c.j])  # contribution to MF
              fnew[c.j] += Av(c.g,mnew[c.i])  # contribution to MF
          return fnew
        fold = m2mf(mold) # magnetization to mean field
        for c in sc.exchange: c.g *= p # scale the many-body exchange
        sc.set_fields(lambda i: (1.-p)*fold[i]) # Add new magnetic fields
        sc.gs_energy() # perform calculation
        mnew = np.array(sc.get_magnetization()).transpose() # magnetization
        return mnew,sc # return new magnetization
    if mode=="default": # mixing mode
      while True: # infinite loop
          fnew,sc = get_new(fold) # new magnetization
          fold = mix*fnew + (1.0-mix)*fold # new fields
          error = np.sum(np.abs(fold-fnew)) # error
          if error<maxerror: break
          print("Error",error)
    elif mode=="broyden": # mixing mode
        from scipy.optimize import broyden1,broyden2,anderson,linearmixing
        def fopt(fold): # function to optimize
          fnew,sc = get_new(fold) # new magnetization
          return fnew - fold
#        fnew = broyden1(fopt,fold)
        try: fnew = broyden1(fopt,fold,maxiter=20)
        except: return spinchain_meanfield(sc,m0=fold)
#        fnew = anderson(fopt,fold)
        print("Done")
        fnew,sc = get_new(fnew) # new magnetization
#        print(fnew)

    print("Done")
    return sc # return the spin chain object




def Av(A,v):
    """Apply matrix to a vector and return a vector"""
    o = np.matrix(A)*np.matrix(v).T
    return np.array(o).reshape(3)
