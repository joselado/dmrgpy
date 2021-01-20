import numpy as np
from .edtk.edchain import State
from .mps import MPS


def evolve_WF(h,wf,ts=np.linspace(0.,10.,10),dt=1e-3):
    """Time evolve a wavefunction"""
    if not h.is_hermitian(): raise # only for hermitian Hamiltonians
    if type(wf)==State: # ED version
        out = []
        for t in ts:
            wf0 = wf.copy()
            wf1 = wf.MBO.exponential(1j*t*h,wf0)
            out.append(wf1.copy())
        return out
#        return [wf.MBO.exponential(1j*t*h,wf) for t in ts]
    elif type(wf)==MPS: # DMRG version
        from .mpsalgebra import exponential_dmrg
        wf0 = wf.copy() # copy wavefunction
        wfout = [] # empty list
        t0 = 0.0 # initial time
        for i in range(len(ts)):
            t1 = ts[i] # final time
            dt01 = t1-t0 # time difference
            nt = max([int(dt01/dt),2]) # number of time-steps
            wf1 = exponential_dmrg(wf.MBO,h,wf0,dt=1j*dt01,nt0=nt) # evolve
#            wf1 = (1.+1j*dt01*h)*wf0 # taylor expansion
#            wf1 = first_order_exponential(h,wf0,dt01)
            wf0 = wf1.copy() # store
            wfout.append(wf1.copy()) # store
            t0 = t1+0.0 # new old time
        return wfout
    else: raise # unrecognized option


imaginary_exponential = evolve_WF # just a wrapper exponential

def first_order_exponential(h,wf,dt):
    print(dt)
    norm = np.sqrt(wf.dot(wf).real) # norm
    wf1 = (1.+1j*dt*h)*wf #- dt**2*h*(h*wf)/2. + (1j*dt)**3*h*(h*(h*wf))/6.
#    wf1 = (1.+1j*dt*h)*wf - dt**2*h*(h*wf)
    wf1 = wf1.normalize()
#    return wf1
    return norm*wf1


