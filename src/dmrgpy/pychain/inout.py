import numpy as np
from scipy.sparse import csr_matrix

def save_sparse_csr(filename,array):
    np.savez(filename,data = array.data ,indices=array.indices,
             indptr =array.indptr, shape=array.shape )

def load_sparse_csr(filename):
    loader = np.load(filename)
    return csr_matrix((  loader['data'], loader['indices'], loader['indptr']),
                         shape = loader['shape'])



def write(x,y,output_file="DATA.OUT",comment=None):
  fo = open(output_file,"w")
  if comment is not None: fo.write("# "+comment+"\n")
  for (ix,iy) in zip(x,y):
    fo.write(str(ix)+"    "+str(iy)+"\n")
  fo.close()


import time


def write_status(i,indict,outf):
  """Write several things in the output files"""
  # selected energy
  outf["energy"].write(str(i)+"   "+str(indict["energy"])+"\n")
  # all the energies
  outf["energies"].write(str(i)+"    ")
  for e in indict["energies"]: # loop over energies
    outf["energies"].write(str(e)+"   ")
  outf["energies"].write("\n")
  # entropy
  outf["entropy"].write(str(i)+"   "+str(indict["entropy"])+"\n")
  # time
  outf["times"].write(str(i)+"   "+str(time.clock())+"\n")
  # correlation
  outf["correlation"].write(str(i)+"   "+str(indict["correlator"])+"\n")
  # length
  outf["total_length"].write(str(i)+"   "+str(indict["total_length"])+"\n")
  outf["length"].write(str(i)+"   "+str(indict["length"])+"\n")
  for key in outf: outf[key].flush() # write on files



def output_files():
  """Return a dictionary for the output files"""
  outf = dict() # dictionary
  outf["entropy"] = open("ENTROPY.OUT","w") # file for the entropy
  outf["times"] = open("TIME.OUT","w") # file for the entropy
  outf["energy"] = open("ENERGY.OUT","w") # file for the entropy
  outf["energies"] = open("ENERGIES.OUT","w") # file for the entropy
  outf["length"] = open("LENGTH.OUT","w") # file for the entropy
  outf["total_length"] = open("TOTAL_LENGTH.OUT","w") # file for the entropy
  outf["correlation"] = open("CORRELATION.OUT","w") # file for the entropy
  return outf
