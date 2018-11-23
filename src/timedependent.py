from __future__ import print_function
import operatornames
import taskdmrg
import numpy as np

def evolution(self,mode="DMRG",**kwargs):
    if mode=="DMRG":  return evolution_dmrg(self,**kwargs)
    if mode=="ED":  
        from dmrgpy2pychain.timedependent import evolution as evolution_exact
        return evolution_exact(self,**kwargs)

def evolution_dmrg(self,name="XX",i=0,j=0,nt=100,dt=0.01):
    namei,namej = operatornames.recognize(self,name)
    task = {"time_evolution":"true",
            "tevol_site_i":str(i),
            "tevol_site_j":str(j),
            "tevol_nt":str(nt),
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




def correlator(self,window=[-1,10],name="ZZ",es=None,dt=0.1,nt=1000):
    """Compute a certain dynamical correlator"""
    (ts,cs) = evolution(self,dt=dt,nt=nt) # get correlator
    cs = cs*np.exp(1j*self.e0*ts) # factor out the phase
    ss = np.fft.fft(cs) # fourier transform
    ws = np.fft.fftfreq(len(cs),d=dt)*8.*np.pi # fourier frequencies
    from scipy.interpolate import interp1d
    fr = interp1d(ws,ss.real,fill_value=0.0,bounds_error=False)
    fi = interp1d(ws,ss.imag,fill_value=0.0,bounds_error=False)
    if es is None:
        es = np.linspace(window[0],window[1],800)
    gr = fr(es)+ 1j*fi(es) # advanced
    ga = np.conjugate(gr) # retarded
#    gp = fr(es) - fr(-es) + 1j*fi(es) + 1j*fi(-es)
    return (es,gr/np.sqrt(nt))






