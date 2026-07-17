import numpy as np
import scipy.linalg as lg
from .entropytk.correlationentropy import get_correlation_matrix
from .entropytk.correlationentropy import get_four_correlation_tensor


def get_correlation_entropy_density(self,**kwargs):
    """Compute the density of spinless correlation entropy"""
    from . import fermionchain
    cm = get_correlation_matrix(self,**kwargs)
    ce,ws = lg.eigh(cm) # return eigenvalues
    ws = ws.T # transpose
    # per-mode entanglement entropy of a free-fermion/orbital occupation
    # number n: -[n*log(n) + (1-n)*log(1-n)] (Peschel's formula), not just
    # -n*log(n) -- the single-term version silently discarded the
    # complementary contribution, undercounting the entropy by roughly 2x
    # for eigenvalues away from 0/1 (confirmed against the ED backend on
    # examples/entanglement_entropy/correlation_entropy's model). Clipping
    # both ends avoids log(0) without needing the old ce<1e-6 special case,
    # since both terms already vanish smoothly as n->0 or n->1.
    ce = np.clip(ce.real,1e-12,1.-1e-12)
    ss = -(ce*np.log(ce) + (1.-ce)*np.log(1.-ce)) # entropies
    out = 0.*ws[0].real # initialize
    for i in range(len(ss)): out = out + ss[i]*np.abs(ws[i])**2
    if type(self)==fermionchain.Fermionic_Chain: pass
    elif type(self)==fermionchain.Spinful_Fermionic_Chain:
      out = [out[2*i]+out[2*i+1] for i in range(len(out)//2)]
    return np.array(out) # return the density

 

def get_correlation_eigenvalues(self,**kwargs):
    """Return the correlation eigenvalues"""
    cm = get_correlation_matrix(self,**kwargs)
    ce = lg.eigvalsh(cm) # return eigenvalues
    # a sanity check
    if np.max(ce)>1.01: raise # they should be at most 1.
    if np.min(ce)<-0.01: raise # they should be zero the smallest
    return ce


def get_correlation_entropy(self,**kwargs):
    """Return the correlation entropy"""
    vs = get_correlation_eigenvalues(self,**kwargs)
    # see get_correlation_entropy_density's comment: use the symmetric
    # -[v*log(v) + (1-v)*log(1-v)] formula, not just -v*log(v)
    vs = np.clip(vs.real,1e-12,1.-1e-12)
    return -np.sum(vs*np.log(vs) + (1.-vs)*np.log(1.-vs))



def get_correlation_entropy_from_wf(self,**kwargs):
    """Return the correlation entropy of a wavefunction"""
    if self.MBO is not None: # if it has an MBO
        return get_correlation_entropy(self.MBO,wf=self,**kwargs)
    else: raise





def get_correlated_orbitals(self,ordered=True,**kwargs):
    """Return the most correlated orbitals"""
    cm = get_correlation_matrix(self,**kwargs)
    es,vs = lg.eigh(cm) # diagonalize
    vs = np.abs(vs.T)**2 # transpose
    esa = np.abs(es-.5)
    if ordered: vs = [y for (x,y) in sorted(zip(esa,vs),key=lambda x: x[0])]
    return vs


def get_correlated_density(self,n=4,**kwargs):
    """Return the most correlated orbitals"""
    vs = get_correlated_orbitals(self,**kwargs)
    return np.sum(vs[0:n],axis=0)



def get_second_order_correlation_entropy(self,**kwargs):
    """Return the second_order correlation entropy"""
    from .entropytk.correlationentropy import get_highorder_correlation_matrix
    cm = get_highorder_correlation_matrix(self,**kwargs)
#    print(np.sum(np.abs(cm-np.conjugate(cm.T))))
    vs = lg.eigvalsh(cm) # return eigenvalues
    # see get_correlation_entropy_density's comment: use the symmetric
    # -[v*log(v) + (1-v)*log(1-v)] formula, not just -v*log(v)
    vs = np.clip(vs.real,1e-12,1.-1e-12)
    out = -np.sum(vs*np.log(vs) + (1.-vs)*np.log(1.-vs))
#    print(np.round(vs,3),np.round(out,2))
    return out
