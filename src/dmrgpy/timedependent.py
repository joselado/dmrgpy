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
    name[0] = name[0].get_dagger()
    self.execute(lambda: name[0].write(name="dc_multioperator_i.in"))
    self.execute(lambda: name[1].write(name="dc_multioperator_j.in"))
    self.execute( lambda : taskdmrg.write_tasks(self)) # write tasks
    self.execute( lambda : self.run()) # run calculation
    cs = self.get_file("TIME_EVOLUTION.OUT").transpose() # time evolution
    ts = np.array([dt*ii for ii in range(nt)]) # times
    return ts,cs[0].real-1j*cs[1] # return



def evolve_and_measure(self,mode="DMRG",**kwargs):
    """Evolve and measure"""
    if mode=="DMRG": return evolve_and_measure_dmrg(self,**kwargs)
    elif mode=="ED": 
        edobj = self.get_ED_obj() # get the ED object
        h = self.hamiltonian # get the ED object
        return tded.evolve_and_measure(edobj,h,**kwargs)



def evolve_and_measure_dmrg(self,operator=None,nt=1000,
        dt=1e-2,wf=None,**kwargs):
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
    if wf is None: wf = self.wf0 # get ground state
    wf.copy(name="psi_evolve_and_measure.mps") # copy wavefunction
    self.execute(lambda: operator.write(name="time_evolution_multioperator.in"))
    self.execute( lambda : taskdmrg.write_tasks(self)) # write tasks
    self.execute( lambda : self.run()) # run calculation
    cs = self.get_file("TIME_EVOLUTION.OUT").transpose() # time evolution
    ts = np.array([dt*ii for ii in range(int(nt))]) # times
    return ts,cs[0].real-1j*cs[1] # return


def evolution_ABA(self,A=None,B=None,mode="DMRG",wf=None,**kwargs):
    """Apply an operator, evolve and measure"""
    if mode=="DMRG":
        if wf is None: wf = self.get_gs() # get ground state
        from .applyoperator import applyoperator
        wfA = applyoperator(self,A,wf) # apply wavefunction 
        return evolve_and_measure_dmrg(self,wf=wfA,operator=B,**kwargs)
    elif mode=="ED":
        edobj = self.get_ED_obj() # get the ED object
        return tded.evolution_ABA(edobj,h=self.hamiltonian,A=A,B=B,wf=wf,
                **kwargs)






def dynamical_correlator(self,window=[-1,10],es=None,dt=0.1,
        nt=None,factor=1,delta=1e-1,**kwargs):
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






