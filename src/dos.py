import numpy as np

def get_dos(mbc,i=0,n=20,delta=0.1):
    """Compute the density of states"""
    mbc.to_folder() # go to the folder
    mbc.setup_sweep("fast")
    mbc.write_hamiltonian() # write the hamiltonian
    task = {"nexcited":str(n),"dos_site":str(i)}
    mbc.setup_task("dos",task=task)
    mbc.run() # run calculation
    es = np.genfromtxt("EXCITED.OUT") # get energies
    es = [es[i]-es[0] for i in range(0,len(es))] # excitation energies
    cs = np.genfromtxt("DOS_ELEMENT.OUT").transpose() # matrix elements
    mbc.to_origin() # go back
    cs = cs[0]**2 + cs[1]**2 # moduli square
    # now convert the cs into a matrix
    cs = cs.reshape((len(es),len(es))) # matrix with the overlaps
    energies = np.linspace(-3.0,3.0,1000) # energies
    dos = np.zeros(energies.shape) + 0j # initizlize
    for i in range(len(dos)): # loop
        for ii in range(len(es)): # loop over states
            for jj in range(len(es)): # loop over states
              dos[i] += cs[ii,jj]/(energies[i] - es[ii] -es[jj] +1j*delta) 
    dos = -dos.imag # imaginary part
    np.savetxt("DOS.OUT",np.matrix([energies,dos]).T) # save
    return (energies,dos)
