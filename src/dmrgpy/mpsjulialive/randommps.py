import numpy as np

from . import taskdmrg
from .mps import MPS

def random_mps_dummy(self,normalize=True):
    """Generate a random MPS, beware of a complicated statistial
    distribution"""
    task = {"random_mps":"true",
            }
    self.task = task
    self.execute( lambda : taskdmrg.write_tasks(self)) # write tasks
    self.execute( lambda : self.run()) # run calculation
    out = MPS(self,name="random.mps").copy() # copy
    norm = np.sqrt(out.overlap(out))
    return (1./norm)*out # return the eigenvector 






def random_product_state(self):
    """Generate a random product state, this ensures a correct
    statistical distribution of the generated states"""
    raise # this is not functional now
    from .fermionchain import Fermionic_Chain, Spinful_Fermionic_Chain
    from .spinchain import Spin_Chain
    mbc = self.clone() # clone the object
    mbc.maxm = 5
    mbc.nsweeps = 10
    h = 0 
    if type(mbc)==Fermionic_Chain or type(mbc)==Spinful_Fermionic_Chain: 
        for i in range(len(mbc.N)): # loop over sites
            if np.random.randint(2)==0: h = h + mbc.N[i] # empty
            else: h = h - mbc.N[i] # full
    elif type(mbc)==Spin_Chain: # spin chain
        for i in range(len(mbc.Sz)): # loop over sites
            m = np.random.random(3) -.5 # random magnetic field
            h = h + m[0]*mbc.Sx[i] + m[1]*mbc.Sy[i] + m[2]*mbc.Sz[i]  # empty
    else: raise # not implemented
    mbc.set_hamiltonian(h) # set the Hamiltonian
    wf = mbc.get_gs() # get the ground state
    wf.set_MBO(self) # set the correct MBO
    return wf


def random_mps(self):
    """Return a random MPS"""
#    try: return random_product_state(self)
#    except:
#        print("Using the default routine")
    wfr = random_mps_dummy(self)
    wfi = random_mps_dummy(self)
    wf = wfr + 1j*wfi # just a linear combination
    wf = wf.normalize()
    return wf


def orthogonal_random_mps(self,wfs):
    """Generate an MPS that is orthogonal to a list"""
    from .algebra.arnolditk import gram_smith_single
    while True: # infinite loop
        wf = self.random_mps() # generate a random MPS
        wf = gram_smith_single(wf,wfs)
        if wf is not None: return wf # return this wavefunction




def random_state(self,mode="DMRG",orthogonal=None,ntries=-1):
    """Generate a random state"""
    if orthogonal is None:
        if mode=="DMRG": return self.random_mps()
        elif mode=="ED": return self.get_ED_obj().random_state()
    else:
        from .algebra.krylov import gram_smith_single
        ii = 0 # counter
        while True: # infinite loop
            wf = random_state(self,mode=mode) # generate a random state
            wf = gram_smith_single(wf,orthogonal)
            if wf is not None: return wf # return this wavefunction
            ii += 1 # increase counter
            if ii==ntries: return None # return anyway




def random_states(self,n=10,orthogonal=True,**kwargs):
    """Return several random states"""
    wfs = []
    for i in range(n):
        if orthogonal:
            wf = random_state(self,orthogonal=wfs,
                    ntries=1,**kwargs) # get state
            if wf is None: break
            else: wfs.append(wf)
        else:
            wfs.append(random_state(self,**kwargs))
    return wfs









