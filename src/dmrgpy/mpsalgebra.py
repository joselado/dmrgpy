from . import mps
import numpy as np

def _use_cpp_ext(self,*wfs):
    """Whether the in-process extension session should be used for a call
    involving the given wavefunction(s). Checking the session alone is not
    enough: not every MPS-producing code path is ported yet (e.g.
    randommps.py's random_mps_dummy() still goes through the old
    file-based backend even when a session is active, since is_hermitian()
    -- called on every gs_energy() -- uses it), so a wavefunction reaching
    here may be file-based (cpp_handle is None) even on a chain with an
    active session. Falling through to the old path in that case is what
    keeps the two backends safely mixable on the same chain."""
    if not (getattr(self,"use_cpp_extension",False) and self._session is not None):
        return False
    return all(wf.cpp_handle is not None for wf in wfs)

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
    """Compute the exponential of a wavefunction"""
    if not self.is_hermitian(h): raise
    if nt0 is None: nt0 = int(h.get_bandwidth(self)*nt)
    # the extension only implements the custom_exp (2nd-order Taylor)
    # variant, which is what tevol_custom_exp=True (the default) selects
    # below in the old-backend path too; fall through to that path if it's
    # been explicitly turned off (toExpH<ITensor> is not ported)
    if self.tevol_custom_exp and _use_cpp_ext(self,wfa):
        tau = complex(-dt.real,dt.imag)
        handle = self._session.exponential_apply(h.to_terms(),wfa.cpp_handle,
                tau,int(nt0))
        return mps.MPS(self,cpp_handle=handle).copy()
    task = {"exponential_eMwf":"true",
            "tevol_dt_real":str(-dt.real),
            "tevol_dt_imag":str(dt.imag),
            "tevol_n":str(int(nt0)),
            }
    if self.tevol_custom_exp: task["tevol_custom_exp"] = "true"
    self.task = task # override tasks
    self.execute(lambda: wfa.write(name="input_wavefunction.mps")) # copy WF
    self.execute(lambda: h.write(name="hamiltonian.in"))
    self.execute(lambda: self.run()) # run calculation
    wf = mps.MPS(self,name="output_wavefunction.mps").copy() # output
    return wf

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
    if _use_cpp_ext(self,wf1,wf2):
        return self._session.overlap(wf1.cpp_handle,wf2.cpp_handle)
    task = {"overlap":"true",
            }
    self.task = task # override tasks
    wf1.write(name="overlap_wf1.mps") # copy wavefunction
    wf2.write(name="overlap_wf2.mps") # copy wavefunction
    self.execute( lambda : self.run()) # run calculation
    m = self.execute( lambda : np.genfromtxt("OVERLAP.OUT")) # run calculation
    return m[0] + 1j*m[1]


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
    if _use_cpp_ext(self,wf1,wf2):
        return self._session.overlap_aMb(wf1.cpp_handle,A.to_terms(),wf2.cpp_handle)
    task = {"overlap_aMb":"true",
            }
    self.task = task # override tasks
    self.execute(lambda: wf1.write(name="overlap_aMb_wf1.mps"))
    self.execute(lambda: wf2.write(name="overlap_aMb_wf2.mps"))
    self.execute(lambda: A.write(name="overlap_aMb_M.in"))
    self.execute(lambda: self.run()) # run calculation
    m = self.execute(lambda : np.genfromtxt("OVERLAP_aMb.OUT")) # read
    return m[0] + 1j*m[1]


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
    if _use_cpp_ext(self,wf1,wf2):
        handle = self._session.sum_mps(wf1.cpp_handle,wf2.cpp_handle)
        return mps.MPS(self,cpp_handle=handle).copy()
    self.execute(lambda: wf1.write(name="summps_wf1.mps")) # write WF
    self.execute(lambda: wf2.write(name="summps_wf2.mps")) # write WF
    task = {"summps":"true",
            }
    self.task = task
    self.execute(lambda : self.run()) # run calculation
    return mps.MPS(self,name="summps_wf3.mps").copy() # copy



def applyoperator_dmrg(self,A,wf):
    """Apply operator to a many body wavefunction"""
    if _use_cpp_ext(self,wf):
        return applyoperator_dmrg_cpp_ext(self,A,wf)
    self.execute(lambda: wf.write(name="applyoperator_wf0.mps")) # write WF
    task = {"applyoperator":"true",
            "applyoperator_wf0":"applyoperator_wf0.mps",
            "applyoperator_multioperator":"applyoperator_multioperator.in",
            "applyoperator_wf1":"applyoperator_wf1.mps",
            }
    self.execute(lambda: A.write(name="applyoperator_multioperator.in"))
    self.task = task
    self.execute( lambda : self.run()) # run calculation
    return mps.MPS(self,name="applyoperator_wf1.mps").copy() # copy


def applyoperator_dmrg_cpp_ext(self,A,wf):
    """Apply operator via the in-process pybind11 extension
    (mpscpp2/chain_session.h's Chain::apply_operator), mirroring
    applyoperator_dmrg() exactly but with no file I/O."""
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    handle = self._session.apply_operator(A.to_terms(),wf.cpp_handle)
    return mps.MPS(self,cpp_handle=handle).copy()


def applyinverse_dmrg(self,A,wf,delta=None,maxn=None):
    """Apply operator to a many body wavefunction"""
    from .algebra.inverse import solve_Ab
    if delta is None: delta = self.cvm_tol # overwrite
    if maxn is None: maxn = self.cvm_nit # overwrite
    if _use_cpp_ext(self,wf):
        self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
        self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
        handle = self._session.apply_inverse(A.to_terms(),wf.cpp_handle,
                delta,int(maxn))
        return mps.MPS(self,cpp_handle=handle).copy()
#    return solve_Ab(A,wf,tol=delta,nmax=1e2)
    self.execute(lambda: wf.write(name="apply_inverse_wf0.mps")) # write WF
    task = {"apply_inverse":"true",
            "cvm_tol":delta,
            "cvm_nit":maxn,
            }
    self.execute(lambda: A.write(name="apply_inverse_multioperator.in"))
    self.task = task
    self.execute( lambda : self.run()) # run calculation
    return mps.MPS(self,name="apply_inverse_wf1.mps").copy() # copy



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
        if self.itensor_version==2:
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
    if _use_cpp_ext(self,wf):
        handle = self._session.conjugate(wf.cpp_handle)
        return mps.MPS(self,cpp_handle=handle).copy()
    self.execute(lambda: wf.write(name="wf1.mps")) # write WF
    task = {"conjugate_mps":"true",
            }
    self.task = task
    self.execute( lambda : self.run()) # run calculation
    return mps.MPS(self,name="wf2.mps").copy() # copy


from .mpsalgebratk.trace import trace
from .mpsalgebratk.trace import inverse_trace


from .mpsalgebratk.disentangle import disentangle_manifold




