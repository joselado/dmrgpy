import numpy as np


def get_dos(mbc,i=0,n=20,delta=0.1,window=5.):
    """Compute the DOS using the KPM"""
    (x1,y1) = mbc.get_dynamical_correlator(i=i,j=i,name="cdc",
            delta=delta,window=[0.0001,window])
    (x2,y2) = mbc.get_dynamical_correlator(i=i,j=i,name="ccd",
            delta=delta,window=[0.0001,window])
    x = np.concatenate([x1,-x2]) # concatenate
    y = np.concatenate([y1,y2]) # concatenate
#    y = [iy for (ix,iy) in sorted(zip(x,y))]
    from scipy.interpolate import interp1d
    f = interp1d(x, y.real,fill_value=0.0,bounds_error=False)
    ne = int(window/delta)*4 # number of energies
    x = np.linspace(-window,window,ne) # new array
    y = f(x) # interpolate

#    from scipy.signal import savgol_filter
#    nw = ne//50
#    if nw%2==0: nw += 1 # odd dnumber
#    y = savgol_filter(y,nw,2) # smooth the curve
    y = convolve(x,y,delta=delta*4)
    np.savetxt("DOS.OUT",np.matrix([x,y]))
    return (x,y)


def convolve(x,y,delta=None):
  """Add a broadening to a DOS"""
  if delta is None: return y # do nothing
  delta = np.abs(delta) # absolute value
  xnew = np.linspace(-1.,1.,len(y)) # array
  d2 = delta/(np.max(x) - np.min(x)) # effective broadening
  fconv = d2/(xnew**2 + d2**2) # convolving function
  yout = np.convolve(y,fconv,mode="same") # same size
  # ensure the normaliation is the same
  ratio = np.sum(np.abs(y))/np.sum(np.abs(yout))
  yout *= ratio
  return yout # return new array




def naive_get_dos(mbc,i=0,n=20,delta=0.1):
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
