import numpy as np
import os




def multicorrelator(self,rs=None,es=np.linspace(0.0,4.0,100),**kwargs):
  """Compute the dynamical corraltor at different energies"""
  if rs is None:
      rs = [[i,0] for i in range(self.ns)] # initialize
  rs = np.array(rs) # to array and transpose
  ds = [] # empty list with total DOS
  for ii in range(self.ns): # loop over spins
    print("Computing",ii)
    (ei,di) = self.get_dynamical_correlator(i=ii,j=ii,es=es,**kwargs)
    ds.append(di) # store
  ds = np.array(ds) 
  os.system("rm -rf MULTILDOS") # remove folder
  os.system("mkdir MULTILDOS") # create folder
  fo = open("MULTILDOS/MULTILDOS.TXT","w") # files with the names
  ie = 0 # start
  for e in es: # loop over energies
    die = np.abs(ds[:,ie]) # correlator at this energy
    ie += 1 # increase
    name0 = "LDOS_"+str(e)+"_.OUT" # name of the output
    np.savetxt("MULTILDOS/"+name0,np.matrix([rs[:,0],rs[:,1],die]).T)
    name = "MULTILDOS/" + name0
    fo.write(name0+"\n") # name of the file
    fo.flush() # flush
  fo.close() # close file
  dt = np.zeros(es.shape,dtype=np.complex128) # initialize
  for di in ds: dt += di # add
  np.savetxt("MULTILDOS/DOS.OUT",np.matrix([es,np.abs(dt)]).T) # save the total
  fo = open("CORRELATORMAP.OUT","w")
  for i in range(self.ns): # loop
      for ie in range(len(es)): # loop over energies
        fo.write(str(i)+"  ") # write site
        fo.write(str(es[ie])+"  ") # write site
        fo.write(str(ds[i,ie].real)+"  ") # write site
        fo.write(str(ds[i,ie].imag)+"\n") # write site
  fo.close() # close file



