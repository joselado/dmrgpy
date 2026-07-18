from . import mps
import numpy as np

def exponential(self,h,wf,mode="DMRG",**kwargs):
    """Compute the exponential"""
    mode = wf.mode # mode of the wavefunction
#    if self.mode is not None: mode = self.mode # redefine
    if mode=="DMRG": 
        if h.is_hermitian(): 
            return exponential_dmrg(self,h,wf,dt=1.0,**kwargs)
        elif h.is_antihermitian(): 
            return exponential_dmrg(self,-1j*h,wf,dt=-1j,**kwargs)
        else:
            print("Warning, using 3rd order taylor expansion mode")
            wf1 = h*wf # apply Hamiltonian
            wf2 = h*wf1 # apply Hamiltonian
            return wf + wf1 + wf2/2.
#            raise
    elif mode=="ED": 
        return self.get_ED_obj().exponential(h,wf,**kwargs)
    else: raise


def exponential_dmrg(self,h,wfa,dt=1.0,nt=1000,nt0=None):
    """Compute the exponential of a wavefunction via the in-process
    pybind11 extension (mpscpp2/chain_session.h's Chain::exponential_apply,
    a custom 2nd-order Taylor expansion)."""
    if not self.is_hermitian(h): raise
    if nt0 is None:
        nt0 = int(h.get_bandwidth(self)*nt)
        # get_bandwidth() runs its own DMRG ground-state search (see
        # bandwidth() in manybodychain.py) and can occasionally
        # underestimate a highly degenerate operator's spectral width
        # depending on the random initial wavefunction; nt0<=0 would divide
        # by a zero/negative effective dt below, so fall back to nt steps
        # (equivalent to a bandwidth of 1) rather than feeding a degenerate
        # step count into the extension.
        if nt0<1: nt0 = nt
    if not self.tevol_custom_exp:
        raise NotImplementedError(
                "tevol_custom_exp=False selects ITensor's toExpH variant, "
                "which only ever existed in the removed file-based backend; "
                "the in-process extension only implements the custom_exp "
                "(2nd-order Taylor) variant, so leave tevol_custom_exp=True")
    tau = complex(-dt.real,dt.imag)
    handle = self._session.exponential_apply(h.to_terms(),wfa.cpp_handle,
            tau,int(nt0))
    return mps.MPS(self,cpp_handle=handle).copy()

def overlap(self,wf1,wf2,mode="DMRG"):
    if self.mode is not None: mode = self.mode # redefine
    if mode=="DMRG": return overlap_dmrg(self,wf1,wf2)
    elif mode=="ED": return self.get_ED_obj().overlap(wf1,wf2)
    else: raise


def overlap_aMb(self,wf1,A,wf2,mode="DMRG"):
    """Compute the overlap <wf1|M|wf2>"""
    if self.mode is not None: mode = self.mode # redefine
    #return wf1.dot(A*wf2) # workaround
    if mode=="DMRG": return overlap_aMb_dmrg(self,wf1,A,wf2)
    elif mode=="ED": return wf1.dot(A*wf2) # workaround
    else: raise


def overlap_dmrg(self,wf1,wf2):
    """Compute the overlap between wavefunctions"""
    return self._session.overlap(wf1.cpp_handle,wf2.cpp_handle)


def overlap_aMb_dmrg(self,wf1,A,wf2):
    """Compute the overlap between wavefunctions"""
    from .multioperator import MultiOperator
    from .multioperatortk.staticoperator import StaticOperator
    if type(A)==StaticOperator:
        return A.aMb(wf1,wf2)
#        return wf1.dot(A*wf2) # workaround
    else:
        return overlap_aMb_dmrg_MO(self,wf1,A,wf2)



def overlap_aMb_dmrg_MO(self,wf1,A,wf2):
    """Compute the overlap between wavefunctions, with A a multioperator"""
    from .multioperator import obj2MO
    A = obj2MO(A) # convert to a MO
    return self._session.overlap_aMb(wf1.cpp_handle,A.to_terms(),wf2.cpp_handle)


def applyoperator(self,A,wf,**kwargs):
    if type(wf)==mps.MPS: mode="DMRG"
    elif type(wf)==np.ndarray: mode="ED"
    else: raise
    if mode=="DMRG": return applyoperator_dmrg(self,A,wf)
    elif mode=="ED": 
        return self.get_ED_obj().applyoperator(A,wf)


def applyinverse(self,A,wf,**kwargs):
    from .edtk.edchain import State
    if type(wf)==mps.MPS: mode="DMRG"
    elif type(wf)==State: mode="ED"
    else: raise
    if mode=="DMRG": return applyinverse_dmrg(self,A,wf,**kwargs)
    elif mode=="ED": 
        return wf.applyinverse(A)
#        return self.get_ED_obj().applyoperator(A,wf)


def summps(self,wf1,wf2,**kwargs):
    if type(wf1)==mps.MPS: mode="DMRG"
    elif type(wf1)==np.ndarray: mode="ED"
    else: raise
    if mode=="DMRG": return summps_dmrg(self,wf1,wf2)
    elif mode=="ED": return wf1 + wf2 #self.get_ED_obj().summps(A,wf1,wf2)



def summps_dmrg(self,wf1,wf2):
    """Apply operator to a many body wavefunction"""
    handle = self._session.sum_mps(wf1.cpp_handle,wf2.cpp_handle)
    return mps.MPS(self,cpp_handle=handle).copy()


def applyoperator_dmrg(self,A,wf):
    """Apply operator via the in-process pybind11 extension
    (mpscpp2/chain_session.h's Chain::apply_operator)."""
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_verbose(self.verbose)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    handle = self._session.apply_operator(A.to_terms(),wf.cpp_handle)
    return mps.MPS(self,cpp_handle=handle).copy()


def applyinverse_dmrg(self,A,wf,delta=None,maxn=None):
    """Apply operator to a many body wavefunction"""
    if delta is None: delta = self.cvm_tol # overwrite
    if maxn is None: maxn = self.cvm_nit # overwrite
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_verbose(self.verbose)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    handle = self._session.apply_inverse(A.to_terms(),wf.cpp_handle,
            delta,int(maxn))
    return mps.MPS(self,cpp_handle=handle).copy()



def operator_norm(self,op,ntries=5,simplify=True):
    """Given a certain operator, compute its norm"""
    if simplify: op = op.simplify() # simplify the operator
    out = [] # empty list
    for i in range(ntries):
        wf = self.random_mps() # random wavefunction
        wf = op*wf # apply the operator
        o = (wf.overlap(wf)).real
        out.append(o)
    return np.mean(out) # return the norm



def is_hermitian(self,op):
    """Given a certain operator, check if it is Hermitian"""
    op = op - op.get_dagger()
    wf = self.random_mps() # random wavefunction
    wf = op*wf # apply the operator
    norm = (wf.dot(wf)).real
    return not norm>1e-4






from .algebra.arnolditk import mpsarnoldi
from .algebra.arnolditk import lowest_energy as lowest_energy_arnoldi
from .algebra.arnolditk import lowest_energy_non_hermitian as lowest_energy_non_hermitian_arnoldi
from .algebra.arnolditk import gram_smith_single


def toMPO(self,H,mode="DMRG"):
    """Transport an operator into a matrix-product operator"""
    if mode=="DMRG":
        if self.itensor_version in (2,3,"python"):
            from .multioperatortk.staticoperator import StaticOperator
            return StaticOperator(H,self) 
        elif self.itensor_version=="julia_live":
            from .mpsjulialive.mpo import MPO
            return MPO(H,MBO=self)
        else: raise # not implemented
    elif mode=="ED":
        from .edtk.edchain import EDOperator
        return EDOperator(H,self.get_ED_obj())
    else: raise



def conjugate_mps(self,wf):
    """Apply operator to a many body wavefunction"""
    handle = self._session.conjugate(wf.cpp_handle)
    return mps.MPS(self,cpp_handle=handle).copy()


from .mpsalgebratk.trace import trace
from .mpsalgebratk.trace import inverse_trace


from .mpsalgebratk.disentangle import disentangle_manifold




