from . import mps
import numpy as np

def exponential(self,h,wf,mode="DMRG",**kwargs):
    """Compute the exponential"""
    if mode=="DMRG": 
        if h.is_hermitian(): 
            return exponential_dmrg(self,h,wf,dt=1.0,**kwargs)
        elif h.is_antihermitian(): 
            return exponential_dmrg(self,-1j*h,wf,dt=1j,**kwargs)
        else:
            print("Operator is not Hermitian nor anti-Hermitian")
            raise
    elif mode=="ED": 
        return self.get_ED_obj().exponential(h,wf,**kwargs)
    else: raise


def exponential_dmrg(self,h,wfa,dt=1.0,nt=1000):
    """Compute the exponential of a wavefunction"""
    if not h.is_hermitian(): raise
    nt0 = int(h.get_bandwidth(self)*nt)
    task = {"exponential_eMwf":"true",
            "tevol_dt_real":str(-dt.real),
            "tevol_dt_imag":str(dt.imag),
            "tevol_n":str(nt0),
            }
    self.task = task # override tasks
    self.execute(lambda: wfa.write(name="input_wavefunction.mps")) # copy WF
    self.execute(lambda: h.write(name="hamiltonian.in"))
    self.execute(lambda: self.run()) # run calculation
    wf = mps.MPS(self,name="output_wavefunction.mps").copy() # output
    return wf

def overlap(self,wf1,wf2,mode="DMRG"):
    if mode=="DMRG": return overlap_dmrg(self,wf1,wf2)
    elif mode=="ED": return self.get_ED_obj().overlap(wf1,wf2)
    else: raise


def overlap_dmrg(self,wf1,wf2):
    """Compute the overlap between wavefunctions"""

    task = {"overlap":"true",
            }
    self.task = task # override tasks
    wf1.write(name="overlap_wf1.mps") # copy wavefunction
    wf2.write(name="overlap_wf2.mps") # copy wavefunction
    self.execute( lambda : self.run()) # run calculation
    m = self.execute( lambda : np.genfromtxt("OVERLAP.OUT")) # run calculation
    return m[0] + 1j*m[1]


def applyoperator(self,A,wf,**kwargs):
    if type(wf)==mps.MPS: mode="DMRG"
    elif type(wf)==np.ndarray: mode="ED"
    else: raise
    if mode=="DMRG": return applyoperator_dmrg(self,A,wf)
    elif mode=="ED": 
        return self.get_ED_obj().applyoperator(A,wf)


def applyinverse(self,A,wf,**kwargs):
    if type(wf)==mps.MPS: mode="DMRG"
    elif type(wf)==np.ndarray: mode="ED"
    else: raise
    if mode=="DMRG": return applyinverse_dmrg(self,A,wf,**kwargs)
    elif mode=="ED": raise
#        return self.get_ED_obj().applyoperator(A,wf)


def summps(self,wf1,wf2,**kwargs):
    if type(wf1)==mps.MPS: mode="DMRG"
    elif type(wf1)==np.ndarray: mode="ED"
    else: raise
    if mode=="DMRG": return summps_dmrg(self,wf1,wf2)
    elif mode=="ED": return wf1 + wf2 #self.get_ED_obj().summps(A,wf1,wf2)



def summps_dmrg(self,wf1,wf2):
    """Apply operator to a many body wavefunction"""
    self.execute(lambda: wf1.write(name="summps_wf1.mps")) # write WF
    self.execute(lambda: wf2.write(name="summps_wf2.mps")) # write WF
    task = {"summps":"true",
            }
    self.task = task
    self.execute(lambda : self.run()) # run calculation
    return mps.MPS(self,name="summps_wf3.mps").copy() # copy



def applyoperator_dmrg(self,A,wf):
    """Apply operator to a many body wavefunction"""
    self.execute(lambda: wf.write()) # write WF
    task = {"applyoperator":"true",
            "applyoperator_wf0":wf.name,
            "applyoperator_multioperator":"applyoperator_multioperator.in",
            "applyoperator_wf1":"applyoperator_wf1.mps",
            }
    self.execute(lambda: A.write(name="applyoperator_multioperator.in"))
    self.task = task
    self.execute( lambda : self.run()) # run calculation
    return mps.MPS(self,name="applyoperator_wf1.mps").copy() # copy


def applyinverse_dmrg(self,A,wf,tol=1e-4,maxn=100):
    """Apply operator to a many body wavefunction"""
    self.execute(lambda: wf.write(name="apply_inverse_wf0.mps")) # write WF
    task = {"apply_inverse":"true",
            "cvm_tol":tol,
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



from .algebra.arnolditk import mpsarnoldi
from .algebra.arnolditk import gram_smith_single
