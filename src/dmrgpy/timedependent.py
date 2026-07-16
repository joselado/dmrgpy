from __future__ import print_function
from . import operatornames
from . import taskdmrg
import numpy as np
import os
from scipy.interpolate import interp1d
from . import multioperator
from .edtk import timedependent as tded



def evolution_DC(self,mode="DMRG",**kwargs):
    if mode=="DMRG":  return evolution_dmrg_DC(self,**kwargs)
    if mode=="ED": 
        edobj = self.get_ED_obj() # get the ED object
        return tded.evolution_DC(edobj,h=self.hamiltonian,**kwargs)



def evolution_dmrg_DC(self,name="XX",nt=10000,dt=0.1,restart=True,**kwargs):
    name = operatornames.str2MO(self,name,**kwargs)
    name[0] = name[0].get_dagger()
    A,B = name[0],name[1]
    if getattr(self,"use_cpp_extension",False) and self._session is not None:
        return evolution_dmrg_DC_cpp_ext(self,A,B,nt,dt)
    if self.fit_td: fittd = "true"
    else: fittd = "false"
    if self.tevol_custom_exp: tevol_custom_exp = "true"
    else: tevol_custom_exp = "false"
    fittd = "true"
    task = {"time_evolution":"true",
            "tevol_nt":str(nt),
            "tevol_fit":fittd,
            "tevol_custom_exp":tevol_custom_exp,
            "tevol_dt":str(dt),
            }
    self.task = task # override tasks
    if restart: # restart the calculation
      self.execute(lambda: os.system("cp psi_GS.mps psi_time_evolution.mps"))
    self.execute(lambda: A.write(name="dc_multioperator_i.in"))
    self.execute(lambda: B.write(name="dc_multioperator_j.in"))
    self.execute( lambda : taskdmrg.write_tasks(self)) # write tasks
    self.execute( lambda : self.run()) # run calculation
    cs = self.get_file("TIME_EVOLUTION.OUT").transpose() # time evolution
    ts = np.array([dt*ii for ii in range(nt)]) # times
    return ts,cs[0].real-1j*cs[1] # return


def evolution_dmrg_DC_cpp_ext(self,A,B,nt,dt):
    """
    Real-time quench dynamical correlator via the in-process pybind11
    extension (mpscpp2/chain_session.h's Chain::quench), mirroring
    evolution_dmrg_DC()'s DMRG path exactly but with no file I/O.

    fit_td is hardcoded False here, not read from self.fit_td: the old
    path writes it to tasks.in under the key "tevol_fit", but
    time_evolution.h actually reads "tevol_fit_td" (a pre-existing,
    unrelated key-name mismatch) -- so the fitApplyMPO branch there is
    unreachable today regardless of self.fit_td, and False reproduces
    that actual behavior rather than the intended-but-never-taken one.

    "restart" is not passed through either: quench()'s C++ implementation
    always starts from get_gs() regardless of its value -- the
    restart/"cp psi_GS.mps psi_time_evolution.mps" step above is
    vestigial in the old backend too, since psi_time_evolution.mps is
    only ever written, never read, by quench().
    """
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    correlator,_wf = self._session.quench(
            self.hamiltonian.to_terms(),A.to_terms(),B.to_terms(),
            int(nt),dt,False)
    cs = np.array(correlator)
    ts = np.array([dt*ii for ii in range(nt)])
    return ts,cs.real-1j*cs.imag



def evolve_and_measure(self,mode="DMRG",**kwargs):
    """Evolve and measure"""
    if mode=="DMRG": return evolve_and_measure_dmrg(self,**kwargs)
    elif mode=="ED": 
        edobj = self.get_ED_obj() # get the ED object
        h = self.hamiltonian # get the ED object
        return tded.evolve_and_measure(edobj,h,**kwargs)



def evolve_and_measure_dmrg(self,operator=None,nt=1000,h=None,
        dt=1e-2,wf=None,**kwargs):
    if h is None: h = self.hamiltonian # Hamiltonian
    if wf is None: wf = self.wf0 # get ground state
    if (getattr(self,"use_cpp_extension",False) and self._session is not None
            and wf.cpp_handle is not None):
        return evolve_and_measure_dmrg_cpp_ext(self,operator,h,wf,nt,dt)
    if self.fit_td: fittd = "true"
    else: fittd = "false"
    if self.tevol_custom_exp: tevol_custom_exp = "true"
    else: tevol_custom_exp = "false"
    fittd = "true"
    task = {"evolution_measure":"true",
            "tevol_nt":str(int(nt)),
            "tevol_fit":fittd,
            "tevol_custom_exp":tevol_custom_exp,
            "tevol_dt":str(dt),
            }
    self.task = task # override tasks
    wf.write(name="psi_evolve_and_measure.mps") # copy wavefunction
    self.execute(lambda: h.write(name="hamiltonian.in"))
    self.execute(lambda: operator.write(name="time_evolution_multioperator.in"))
    self.execute( lambda : taskdmrg.write_tasks(self)) # write tasks
    self.execute( lambda : self.run()) # run calculation
    cs = self.get_file("TIME_EVOLUTION.OUT").transpose() # time evolution
    ts = np.array([dt*ii for ii in range(int(nt))]) # times
    return ts,cs[0].real-1j*cs[1] # return


def evolve_and_measure_dmrg_cpp_ext(self,operator,h,wf,nt,dt):
    """
    Real-time evolution + measurement via the in-process pybind11
    extension (mpscpp2/chain_session.h's Chain::evolve_and_measure),
    mirroring evolve_and_measure_dmrg()'s DMRG path exactly but with no
    file I/O. fit_td hardcoded False for the same reason as
    evolution_dmrg_DC_cpp_ext (see its docstring): the "tevol_fit"/
    "tevol_fit_td" key-name mismatch means the old backend's fitApplyMPO
    branch is unreachable regardless of self.fit_td.
    """
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    correlator,_wf = self._session.evolve_and_measure(
            h.to_terms(),operator.to_terms(),wf.cpp_handle,
            int(nt),dt,False)
    cs = np.array(correlator)
    ts = np.array([dt*ii for ii in range(int(nt))])
    return ts,cs.real-1j*cs.imag


def evolution_ABA(self,A=None,B=None,mode="DMRG",wf=None,**kwargs):
    """Apply an operator, evolve and measure"""
    if A is None: A = multioperator.identity()
    if B is None: B = multioperator.identity()
    if mode=="DMRG":
        if wf is None: wf = self.get_gs() # get ground state
        wfA = A*wf # apply the operator
        return evolve_and_measure_dmrg(self,wf=wfA,operator=B,**kwargs)
    elif mode=="ED":
        edobj = self.get_ED_obj() # get the ED object
        return tded.evolution_ABA(edobj,h=self.hamiltonian,A=A,B=B,wf=wf,
                **kwargs)






def dynamical_correlator(self,window=[-1,10],es=None,dt=0.1,
        nt=None,factor=1,delta=5e-2,**kwargs):
    """Compute a certain dynamical correlator"""
    self.get_gs() # get the ground state
    if nt is None: nt=int(100/delta/dt)
    (ts,cs) = evolution_DC(self,dt=dt,nt=nt,**kwargs) # get correlator
#    cs = cs*np.exp(-1j*self.e0*ts) # factor out the phase
    # interpolate the time evolution
    ftr = interp1d(ts,cs.real,fill_value=0.0,bounds_error=False)
    fti = interp1d(ts,cs.imag,fill_value=0.0,bounds_error=False)
    # interpolate the time evolution
    tnew = np.linspace(np.min(ts),np.max(ts),len(ts)*factor) # ten times
    cnew = ftr(tnew) + 1j*fti(tnew)
    ts = tnew.copy() # overwrite
    cs = cnew.copy() # overwrite
    dtnew = dt/factor
    # do the fourier transform
    ss = np.fft.fft(cs) # fourier transform
    ws = np.fft.fftfreq(len(cs),d=dtnew)*2.*np.pi # fourier frequencies
    fr = interp1d(ws,ss.real,fill_value=0.0,bounds_error=False)
    fi = interp1d(ws,ss.imag,fill_value=0.0,bounds_error=False)
    if es is None:
        es = np.linspace(window[0],window[1],800)
    gr = fr(es)+ 1j*fi(es) # advanced
    ga = np.conjugate(gr) # retarded
#    gp = fr(es) - fr(-es) + 1j*fi(es) + 1j*fi(-es)
    return (es,gr/np.sqrt(nt))




def generic_evolution(H,wf,normalize=True,dt=1e-2,nt=100,A=None):
    """Perform a time evolution and project onto itself,
    assuming U = e^tH """
    wf0 = wf.copy() # copy wavefunction
    wf1 = wf.copy() # copy wavefunction
    out = []
    for i in range(int(nt)): # loop
        wf1 = wf1 + dt*H*wf1
        if normalize:  wf1 = wf1*(1./np.sqrt(wf1.dot(wf1)))
  #      out.append(wf0.dot(wf1)) # compute
        out.append(wf1.dot(A*wf1)) # compute
        print(i)
    return np.array(range(int(nt)))*dt,np.array(out) # retunr result



