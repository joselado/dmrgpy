from __future__ import print_function
import os
import numpy as np
from inout import save_sparse_csr
from inout import load_sparse_csr
from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix


check = True

def read_couplings(cs):
  """Read matrices associated to this couplings"""
  ms = []
  for c in cs:
    name = "coupling_"+str(c[0])+"_"+str(c[1])+"_.op"
    ms.append(read_matrix(name))
  return ms # return matrices


def read_all(path=""):
  """Read all the matrices of the spin hamiltonian"""
  files = os.listdir(path)
  fops = [] # files with operators
  for f in files:  # retain operator files
    if ".op" in f: fops.append(f)
  coupling = dict() # dictionary for the couplings
  zeeman = dict() # dictionary for zeeman
  for f in fops: # loop over operators
    if "coupling_" in f: # coupling matrix
      c = f.split("_")
      i = int(float(c[1])) # index
      j = int(float(c[2])) # index
      coupling[i,j] = read_matrix(f) # read this matrix
    if "zeeman" in f: # zeeman matrix
      c = f.split("_")
      if c[1] == "sx": d = "x"
      if c[1] == "sy": d = "y"
      if c[1] == "sz": d = "z"
      i = int(float(c[2])) # index
      zeeman[d,i] = read_matrix(f) # read this matrix
  class Terms(): pass
  terms = Terms() # create object
  terms.coupling = coupling # store couplings
  terms.zeeman = zeeman # store zeeman
  return terms


def read_everything(path=""):
  """Read all the matrices, and put in a diccionary"""
  files = os.listdir(os.getcwd()+"/"+path)
  fops = []
  for f in files:  # retain operator files
    if ".op" in f: fops.append(f)
  ops = []
  for f in fops:
    ops.append(read_matrix(f))  # append the matrix
  return ops # return list with matrices



def convert_matrix(namefile):
  """ Function which reads a sparse matrix, and saves as .npz"""
  d = open(namefile,"r").readlines()[0] # first line
  d = int(d.split("=")[1]) # dimension of the matrix
  from scipy.sparse import csc_matrix
#  try:
  if True:
#    print "Reading",namefile
    m = np.genfromtxt(namefile) # read the file
    try:
      a = m[0][0] # check if it is a matrix
      m = m.transpose() # transpose matrix, if more than one row
    except:
      m = [[m[0]],[m[1]],[m[2]],[m[3]]] # fix for one row
    row = m[0].astype(int) # row
    col = m[1].astype(int) # column
    data = m[2] + 1j*m[3] # data
    mout = csc_matrix((data,(row,col)),shape=(d,d),dtype=np.complex) # create the matrix
#    if np.max(np.abs(mout.todense() - mout.todense().H))>0.00001: raise
    write_matrix(namefile,mout) # save matrix as numpy object
    return mout

def write_matrix(namefile,m):
  """Write matrix in a file"""
  save_sparse_csr(namefile+".npz",csr_matrix(m)) # save the matrix



def read_matrix(namefile):
  """Read matrix from a file"""
  return csc_matrix(load_sparse_csr(namefile+".npz")) # read the matrix





def convert_sop(spins):
  """Transfrom from C++ format to python format"""
  for i in range(len(spins)):
    convert_matrix("sx_"+str(i)+"_.op") # convert sx
    convert_matrix("sy_"+str(i)+"_.op") # convert sx
    convert_matrix("sz_"+str(i)+"_.op") # convert sx




def read_sop(spins,check=check):
  """ Read all the spin operators"""
  class Sclass: pass
  sc = Sclass() # empty class
  sxi,syi,szi = [],[],[] # empty lists
  for i in range(len(spins)):
    sxi.append(read_matrix("sx_"+str(i)+"_.op")) # read sx
    syi.append(read_matrix("sy_"+str(i)+"_.op")) # read sy
    szi.append(read_matrix("sz_"+str(i)+"_.op")) # read sz
  # now save in the class
  sc.sxi = sxi
  sc.syi = syi
  sc.szi = szi
  # and calculate total spins
  sc.sx = sum(sxi)
  sc.sy = sum(syi)
  sc.sz = sum(szi)
  sc.s2 = sc.sx*sc.sx + sc.sy*sc.sy + sc.sz*sc.sz
  import checking
  for i in range(len(spins)):
    checking.angular(sc.sxi[i],sc.syi[i],sc.szi[i])
    checking.angular(sc.syi[i],sc.szi[i],sc.sxi[i])
    checking.angular(sc.szi[i],sc.sxi[i],sc.syi[i])

  return sc # return spin class






