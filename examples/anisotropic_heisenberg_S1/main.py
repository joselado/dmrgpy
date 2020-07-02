
import numpy as np
from dmrgpy import spinchain
n = 200
spins = ["S=1/2"]+["S=1" for i in range(n)]+["S=1/2"] # spin 1/2 heisenberg chain

def get(D):
    sc = spinchain.Spin_Chain(spins) # create the spin chain
    h = 0
    for i in range(len(spins)-1):
        h = h +sc.Sx[i]*sc.Sx[i+1]
        h = h +sc.Sy[i]*sc.Sy[i+1]
        h = h +sc.Sz[i]*sc.Sz[i+1]
    for i in range(len(spins)):  h = h - D*sc.Sz[i]*sc.Sz[i]
    # define a Neel Hamiltonian
    h0 = 0
    for i in range(len(spins)): h0 = h0 + ((-1)**i)*sc.Sz[i]
    sc.nsweeps = 5
    sc.maxm = 20
    sc.set_hamiltonian(h0) ; wf = sc.get_gs() # get the Neel wavefunction
    
    sc.set_hamiltonian(h)
    nc = len(spins)//3
    wf = sc.get_gs(wf0=wf)
    print(D)
    return wf.get_entropy(len(spins)//2)



xs = np.linspace(0.,.6,30)
ys = [get(d) for d in xs]

np.savetxt("SWEEP.OUT",np.array([xs,ys]).T)

