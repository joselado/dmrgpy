import scipy.sparse.linalg as slg
from ..algebra import kpm
import numpy as np
from scipy.interpolate import interp1d
from .edchain import State

def get_distribution(self,X=None,wf=None,**kwargs):
    """Get a certain distribution"""
    if wf is None: wf = self.get_gs() # ground state
    elif type(wf)==State: wf = wf.v
    X = self.get_operator(X)
    return distribution_kpm(wf,X=X,**kwargs)


def distribution_kpm(wf0,X=None,scale=10.0,
        delta=1e-1,xs=None):
    """Compute <0| \delta (m-M) |0> using the KPM"""
    if X is None: raise
    M = X # assign
    vi = wf0 # first wavefunction
    vj = wf0 # second wavefunction
    emax = slg.eigsh(M,k=1,ncv=20,which="LA")[0] # upper energy
    emin = slg.eigsh(-M,k=1,ncv=20,which="LA")[0] # upper energy
    scale = np.max(np.abs([emin,emax]))*2. # compute the scale
    n = int(scale/delta) # number of polynomials
    (xs2,ys2) = kpm.dm_vivj_energy(M,vi,vj,scale=scale,
                                npol=n*4,ne=n*10)
    ys2 /= np.pi # normalization from the KPM
    if xs is None: return xs2,ys2
    ys = interp1d(xs2,ys2.real,fill_value=0.,bounds_error=False)(xs) 
    ys = ys+ 1j*interp1d(xs2,ys2.real,fill_value=0.,bounds_error=False)(xs) 
    return xs,ys # return correlator



