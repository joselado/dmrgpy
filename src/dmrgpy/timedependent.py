from __future__ import print_function
from . import operatornames
from . import taskdmrg
import numpy as np
from scipy.interpolate import interp1d

def evolution(self,mode="DMRG",**kwargs):
    if mode=="DMRG":  return evolution_dmrg(self,**kwargs)
    if mode=="ED":  
        from .dmrgpy2pychain import timedependent
        evolution_exact = timedependent.evolution
        return evolution_exact(self,**kwargs)

def evolution_dmrg(self,name="XX",i=0,j=0,nt=10000,dt=0.1):
    namei,namej = operatornames.recognize(name)
    namei = operatornames.hermitian(namei) # get the Hermitian one
    if self.fit_td: fittd = "true"
    else: fittd = "false"
    fittd = "true"
    task = {"time_evolution":"true",
            "tevol_site_i":str(i),
            "tevol_site_j":str(j),
            "tevol_nt":str(nt),
            "tevol_fit":fittd,
            "tevol_dt":str(dt),
            "tevol_operator_i":namei,
            "tevol_operator_j":namej,
            }
    self.task = task # override tasks
    self.execute( lambda : taskdmrg.write_tasks(self)) # write tasks
    self.execute( lambda : self.run()) # run calculation
    cs = self.get_file("TIME_EVOLUTION.OUT").transpose() # time evolution
    ts = np.array([dt*ii for ii in range(nt)]) # times
    return ts,cs[0].real-1j*cs[1] # return




def dynamical_correlator(self,window=[-1,10],name="ZZ",es=None,dt=0.1,
        nt=None,factor=1,i=0,j=0,delta=1e-1):
    """Compute a certain dynamical correlator"""
    self.get_gs() # get the ground state
    if nt is None: nt=int(100/delta/dt)
    (ts,cs) = evolution(self,dt=dt,nt=nt,i=i,j=j,name=name) # get correlator
    cs = cs*np.exp(-1j*self.e0*ts) # factor out the phase
    # interpolate the time evolution
    ftr = interp1d(ts,cs.real,fill_value=0.0,bounds_error=False)
    fti = interp1d(ts,cs.imag,fill_value=0.0,bounds_error=False)
    # interpolate the time evolution
#    factor = 10 # interpolation factor
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






