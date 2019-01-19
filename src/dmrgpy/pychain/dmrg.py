from __future__ import print_function
from __future__ import division
from . import tensorialf90 # fortran90 library
from . import traceoverf90
from . import spectrum
from scipy.sparse import linalg as slg # linear algebra library
from scipy import linalg as lg # linear algebra library
import scipy.sparse as sp
from scipy.sparse import csc_matrix
from scipy.sparse import coo_matrix
import numpy as np
from copy import deepcopy
import time
from . import dmrgmethods
from scipy.sparse.linalg import LinearOperator
from . import ppdmrg # posprocessing library
from . import dmrgtk
import os

silent = True

# routines to perform DMRG in spin chains

from .tensorial import tensorial_LO
from . import tensorial


def test_tensorial(n=8):
  r1 = np.random.random((n,n))
  r2 = np.random.random((n,n))
  op1 = tensorial_operator(r1,r2,sparse=False)
  op2 = tensorial_operator(r1,r2,sparse=True).todense()
  if np.max(np.abs(op1-op2))>0.000001: raise

# load needed functions
tensorial_operator = dmrgtk.tensorialoperator
trace_over = dmrgtk.traceover
get_entropy = dmrgtk.getentropy
ground_state = dmrgtk.groundstate
coupled_hamiltonian = dmrgtk.coupledhamiltonian



def states_dmat(dmat,datadict):
  """Most important states in the density matrix"""
  nstates = datadict["number_of_states"] # states to retain
  told = time.clock() # old time
  es,vecs = lg.eigh(-dmat) # call lapack
  es = -es # set in positive
  vecs = np.transpose(vecs) # transpose vectors
  if np.sum(np.abs(dmat-dmat.H))>0.0001: raise # check that it is hermitian
  #  print("Diagonalizing density matrix, dimension",dmat.shape)
  if len(vecs)>nstates: 
    adaptive = datadict["ensure_symmetry"] # cut by value
    if adaptive: # adaptive algorithm
      vouts = [] # empty list
      dmin = es[nstates]*0.99 # minimum density
      for (v,d) in zip(vecs,es): # loop over densities
        if np.abs(d)>np.abs(dmin): # if sizeble contribution
          vouts.append(v)
    else: 
      vouts = [vecs[i] for i in range(nstates)] 
  else: vouts = vecs # not enoght vectors
  if np.min(es)<-0.00000001: print(es) ; raise
  nstates_out = len(vouts) # number of retained states
  if not silent: print("Number of states is DM",nstates_out)
  error = 1. - np.sum(es[0:nstates_out])
  if not silent: print("Truncation error is",error)
  entropy = get_entropy(es) # get the entropy
  if not silent: print("Entropy is",entropy)
  return csc_matrix(np.array((vouts)).T),entropy




def project(op,R):
  """Project an operator"""
  op = csc_matrix(op)
  R = csc_matrix(R)
  return R.H*op*R


def dmrgdict():
  """Returns a dictionary for the dmrg calculation"""
  ddict = dict() # create dictionary
  ddict["target"] = 0 # state(s) to use for the DM
  ddict["diag_states"] = 10 # number of minimum states to diagonalize
  ddict["diag_mode"] = "arpack" # use arpack
  ddict["tol"] = 0.00001 # precission in the diagonalization
  ddict["dynamic_target"] = None # ground state
  ddict["DM_target"] = None # not set
  ddict["ensure_symmetry"] = False # truncate the DM without breaking symmetry
  ddict["retain_states"] = 1 # one state
  ddict["dynamic_retain_states"] = None  # retained states
  ddict["periodic"] = False # open boundary
  ddict["finite"] = False # open boundary
  ddict["number_of_states"] = 30 # retained states
  ddict["dynamic_number_of_states"] = None  # retained states
  ddict["LO"] = False # use linear operator
  ddict["finite_num_ite"] = 2 # number of iterations in finite DMRG
  ddict["avoid_edge"] = 0 # stop before reaching the edge
  ddict["target_length"] = None # stop before reaching the edge
  ddict["target_function"] = None # Function that decides which states to target
  ddict["ndmrg"] = 100
  return ddict


def update_parameters(dmrgp,datadict,idmrg=0):
  """Update the parameters of the DMRG calculation"""
  names = ["number_of_states","retain_states","target"] # different parameters
  for name in names: # loop
    if callable(datadict["dynamic_"+name]): # if it is a function
      dmrgp[name] = datadict["dynamic_"+name](dmrgp,idmrg=idmrg)
    else: # otherwise
      dmrgp[name] = datadict[name]



def onfly_update(dmrgp):
  """Update the dictionary by checking a file on the fly"""
  


class DMRGresult():
    """Class for the result of a DMRG calculation"""
    def __init__(self):
        self.path = os.getcwd()+"/.dmrgfolder/"
        self.iteration_energy = [] # empty list
        self.iteration_energies = [] # empty list
        self.iteration_entropies = [] # empty list
        self.iteration_dmdis = [] # empty list
        os.system("rm -rf "+self.path) # remove temporal folder
        os.system("mkdir "+self.path) # create temporal folder
        self.inipath = os.getcwd() # initial path
    def to_folder(self): os.chdir(self.path)
    def to_origin(self): os.chdir(self.inipath) # go to original folder
    def write(self):
        """Write several results in a file"""
        # write the energy #
        m = self.iteration_energy
        m = np.array([range(len(m)),m]).T
        np.savetxt("ITERATION_ENERGY.OUT",m)
        # write the energies #
        m = self.iteration_energies
        f = open("ITERATION_ENERGIES.OUT","w")
        for i in range(len(m)):
          f.write(str(i)+"   ")
          for j in range(len(m[i])):
              f.write(str(m[i][j])+"  ") # 
          f.write("\n") # 
        f.close()

        # write the DM distance #
        m = self.iteration_dmdis
        f = open("ITERATION_DMDIS.OUT","w")
        for i in range(len(m)):
            f.write(str(i)+"  ") # 
            for j in range(len(m[i])):
              f.write(str(m[i][j])+"  ") # 
            f.write("\n") # 
        f.close()





def infinite_dmrg(datadict):
  """Main DMRG loop"""
  # extract variables
  dmrgout = DMRGresult() # create class
  dmrgout.to_folder() # go to temporal folder
  right = datadict["right"]
  left = datadict["left"]
  onsite = datadict["onsite"]
  periodic = datadict["periodic"]
  finite = datadict["finite"]
  LO = datadict["LO"]
  site_operators = datadict["site_operator_generator"]
  retain_states = datadict["retain_states"]
  ndmrg = datadict["ndmrg"] # number of dmrg steps
  # generate initial guess
  if periodic: raise 
#dmrgp = initial_BoBo(right,left,onsite,site_operators=site_operators)
  else: 
      coupling = datadict["coupling"] # function returning list with couplings
      cs = coupling(0,1) # list with couplings
      dmrgp = initial_BooB(right,left,onsite,
          site_operators=site_operators,cs=cs) 
  for key in datadict: dmrgp[key] = datadict[key] # copy dictionary
  outf = output_files() # get the dictionary
  #########################################
  #########################################
  #########################################
  ######## INFINITE DMRG ALGORITHM ########
  #########################################
  #########################################
  #########################################
  if not finite: # infinite DMRG scheme
    dmrgp["store_operators"] = False # do not store the operators
    for idmrg in range(1,ndmrg): # loop over dmrg steps
      if not silent: 
        print()
        print("##### Iteration number",idmrg,"######")
      dmrgp = fusechain(dmrgp,dmrgp) # create the full chain
      # perform the DMRG step
      dmrgp = dmrg_step(dmrgp,periodic=periodic,LO=LO,dmrgout=dmrgout) 
      update_parameters(dmrgp,datadict,idmrg=idmrg) # update dictionary
      # write all the energies
      write_status(idmrg,dmrgp,outf) # write several things
      # next iteration
#    ppdmrg.correlator(dmrgp) # write results
  #########################################
  #########################################
  #########################################
  ######## FINITE DMRG ALGORITHM ##########
  #########################################
  #########################################
  #########################################
  if finite: # finite DMRG scheme
    idmrg = 0 # start
    length = datadict["target_length"] # length opf the chain
    if length%2 != 0: raise
    ndmrg = length - 4 # number of dmrg steps
    dmrgp["store_operators"] = False # do not store the operators
    nsites2 = ndmrg//2 # number of sites
    dicts_right = [None for i in range(ndmrg)] # empty list to store
    dicts_left = [None for i in range(ndmrg)] # empty list to store
    dicts = [None for i in range(ndmrg)] # empty list to store
    if not silent: print("Initial half-sweep")
    # The length will be ndmrg + 4
    dicts[0] = copydict(dmrgp) # # first one
    dicts_left[0] = copydict(dmrgp) # # first one
    for i in range(0,nsites2): # do the first half sweep
      dmrgp = fusechain(dmrgp,dmrgp,l=None) # create chain for the next iteration
      dmrgp = dmrg_step(dmrgp,dmrgout=dmrgout) # perform the DMRG step
      write_status(idmrg,dmrgp,outf) # write several things
      idmrg += 1 # increase
      dicts[i+1] = copydict(dmrgp) # save this dictionary (has all the matrices)
    # now we have stored the matrices for the different lengths, perform scf
    # keep growing the block to the right (first time)
    jdict = copydict(dmrgp)
    for j in range(0,nsites2-2): # do the second half sweep
      # take the right block of the the kdict iteration
      if not silent: print("#####  Half sweep  ########")
      jdict = fusechain(jdict,dicts[nsites2-j-2]) # create chain for the next iteration
      # now perform the DMRG step
      jdict = dmrg_step(jdict,dmrgout=dmrgout) # perform the DMRG step
      write_status(idmrg,jdict,outf) # write several things
      idmrg += 1 # increase
      dicts[nsites2+j+1] = copydict(jdict) # save this dictionary
    # The first sweep has been completed
    if not silent: 
      print("################################")
      print("First sweep completed")
      print("################################")
    # Do a complete sweep to the left, so that everything gets its
    # right simension
    L = nsites2*2 # total number sweeps
    num_ite = datadict["finite_num_ite"]
    # temporal workaround
    dicts_right = [copydict(i) for i in dicts] # empty list to store
    dicts_left = [copydict(i) for i in dicts] # empty list to store
    for ite in range(num_ite):
      avoid = dmrgp["avoid_edge"]
      if not silent:
        print("################################")
        print("Performing iteration number",ite)
        print("################################")
      ######################
      # start on the right #
      ######################
      jdict = copydict(dicts_right[avoid]) # start over from the beggining
      for i in range(avoid,L-2-avoid): # do a full sweep to the left
#        jdict = fusechain(jdict,dicts[L-i-2])
        if not silent: print("#####  Left sweep  ########")
        jdict = reflectdict(jdict) # reflect the Hamiltonian
        jdict = fusechain(jdict,dicts_left[L-i-2])
        jdict = dmrg_step(jdict,dmrgout=dmrgout) # perform the DMRG step
        write_status(idmrg,jdict,outf) # write several things
        idmrg += 1 # increase
        dicts_right[i+1] = copydict(jdict)
      ######################
      # start on the left #
      ######################
      jdict = copydict(dicts_left[avoid])
      bulkn = (L-2)//4*3 # index for bulk calculation
      if ite==num_ite-1: jdict["store_operators"] = True # store the operators
      jdict["store_full_operators"] = True
      for i in range(avoid,L-2-avoid): # do a full sweep to the right
        if not silent: print("#####  Right sweep  ########")
        jdict = fusechain(jdict,dicts_right[L-i-2])
        jdict = dmrg_step(jdict,dmrgout=dmrgout) # perform the DMRG step
        write_status(idmrg,jdict,outf) # write several things
        idmrg += 1 # increase
        dicts_left[i+1] = copydict(jdict) # store the result
        if ite==num_ite-1: # last iteration
          if i==bulkn-1: # in the next iteration store the full operators
            jdict["store_full_operators"] = True
          if i==bulkn: # when the length is 3/4 of the total
            ppdmrg.correlator(jdict,outfile="CORRELATORS_BULK.OUT") # write results
#            ppdmrg.dynamical_correlator(jdict) # write results
  
    # write results
#  dmrgout.energy = np.genfromtxt("ENERGY.OUT").transpose()[1][-1]
  dmrgout.write() # write all the results
  dmrgout.to_origin() # go to original folder
  for f in outf: 
      outf[f].close()
  return dmrgout

  

def dmrg_step(dmrgp,periodic=False,LO=True,dmrgout=None):
  """Perform a single DMRG step"""
  if periodic: return dmrg_BoBo(dmrgp,LO=LO) # infinite DMRG method
  else: return dmrg_BooB(dmrgp,LO=LO,dmrgout=dmrgout) # infinite DMRG method



from .inout import write_status,output_files


def copydict(dmrgp):
  """Copy elements of a dictionary, at least the possible ones"""
  return  deepcopy(dmrgp)
  outdict = dict()
  try:
    for key in dmrgp:
      try: outdict[key] = deepcopy(dmrgp[key])
      except: 
        if not silent: print("WARNING, not copied",key)
  except: pass
  return outdict



def reflectdict(dmrgp,center=0):
  """Return a dictionary of a chain, but with
  the Hamiltonian reflected from left to right 
  in the site center"""
  odict = copydict(dmrgp) # copy the dictionary
  tmpdict = copydict(dmrgp) # copy the dictionary
  # reflect different functions
  odict["onsiste"] = lambda i: tmpdict["onsite"](center-i) # reflect function
  odict["right"] = lambda i: tmpdict["left"](center-i) # reflect function
  odict["left"] = lambda i: tmpdict["right"](center-i) # reflect function
  return odict


def fusechain(in1,in2,l=None):
  """Function that fuses the relevant term of the two hamiltonians"""
  out = copydict(in1) # copy dictionary
  if l is None: l = in1["length"] # length of the block
#  print("Generating site ",l)
  out["lS"] = in1["onsite"](l) # site
  out["rS"] = in1["onsite"](l+1) # site
  out["rB"] = in2["lB"] # block hamiltonian
  out["rS_rB"] = in2["lS_lB"] # block hamiltonian
  out["rB_rS"] = in2["lB_lS"] # coupling
  out["lS_rS"] = in1["right"](l) # coupling
  out["rS_lS"] = in1["left"](l+1) # coupling
  if in1["store_operators"] and in2["store_operators"]: # if it has operators
    out["left_site_operators"] = in1["site_operator_generator"](l)
    out["right_site_operators"] = in2["site_operator_generator"](l+1)
  if callable(in1["target_function"]): # if there is a function to decide the wave
    out["left_site_operators"] = in1["site_operator_generator"](l)
    out["right_site_operators"] = in2["site_operator_generator"](l+1)
    out["sum_left_block_operators"] = in1["sum_block_operators"]
    out["sum_right_block_operators"] = in2["sum_block_operators"]
  out["total_length"] = in1["length"] + in2["length"] + 2 # total length of the chain
  return out



# different types of DMRG methods


def dmrg_BooB(dmrgp,integrate_right=True,LO=True,dmrgout=None):
    """Perform a single DMRG step,
       geometry is [ ]oo[ ],
       with open boundary conditions"""
    outdict = copydict(dmrgp) # output dictionary
    ##### Couplings #####
    rS_rB = dmrgp["rS_rB"] # coupling in the right site
    lS_lB = dmrgp["lS_lB"] # coupling in the left site
    lS_rS = dmrgp["lS_rS"] # coupling in the left site
    rS_lS = dmrgp["rS_lS"] # coupling in the left site
    rB_rS = dmrgp["rB_rS"] # coupling to the right bloxk
    lB_lS = dmrgp["lB_lS"] # coupling to the left block
    rB = dmrgp["rB"] # right hamiltonian
    lB = dmrgp["lB"] # left hamiltonian
    rS = dmrgp["rS"] # left hamiltonian
    lS = dmrgp["lS"] # left hamiltonian
    coupling = dmrgp["coupling"] # function returning list with couplings
    # perform the calculation
    dim_lS = lS.shape[0] # left site dimension
    dim_rS = rS.shape[0] # right site dimension
    dim_rB = rB.shape[0] # dimension of the block 
    dim_lB = lB.shape[0] # dimension of the block 
    # idensity operators
    id_lB = sp.eye(dim_lB,dtype=np.complex) # identity operator
    id_rB = sp.eye(dim_rB,dtype=np.complex) # identity operator
    id_lS = sp.eye(dim_lS,dtype=np.complex) # identity operator
    id_rS = sp.eye(dim_rS,dtype=np.complex) # identity operator
    id_rSB = sp.eye(dim_rS*dim_rB,dtype=np.complex) # right identity operator
    id_lSB = sp.eye(dim_lS*dim_lB,dtype=np.complex) # right identity operator
    # function to upgrade an operator to the full space
    def m2full(A,mtype=None,iterative=False):
      """Function that builds a full operator"""
      if iterative: return [m2full(Ai,mtype=mtype) for Ai in A]
      if mtype=="lB": # left block
        B = tensorial_operator(id_lS,A)
        B = tensorial_operator(B,id_rSB,LO=LO) 
      elif mtype=="lS": # left block
        B = tensorial_operator(A,id_lB)
        B = tensorial_operator(B,id_rSB,LO=LO) 
      elif mtype=="rB": # left block
        B = tensorial_operator(id_rS,A)
        B = tensorial_operator(id_lSB,B,LO=LO) 
      elif mtype=="rS": # left block
        B = tensorial_operator(A,id_rB)
        B = tensorial_operator(id_lSB,B,LO=LO)
      else: raise
      return B
    # create the full operators if needed
    if callable(dmrgp["target_function"]): # if there is a function
      outdict["left_site_operators_full"] = m2full(dmrgp["left_site_operators"],
                                             iterative=True,mtype="lS")
      outdict["right_site_operators_full"] = m2full(dmrgp["right_site_operators"],
                                             iterative=True,mtype="rS")
      outdict["sum_right_block_operators_full"] = m2full(dmrgp["sum_right_block_operators"],
                                             iterative=True,mtype="rB")
      outdict["sum_left_block_operators_full"] = m2full(dmrgp["sum_left_block_operators"],
                                             iterative=True,mtype="lB")
    # index of the new site
    length = dmrgp["length"]
    ###########################
    # couple site to right block
    ###########################
    cs = coupling(length+2,length+3) # list with couplings
    rS_plus_rB = (tensorial_operator(id_rS,rB) +
                   tensorial_operator(rS,id_rB) + 
                   coupled_hamiltonian(rS_rB,rB_rS,cs=cs) ) # block+right site 
    ###########################
    # couple site to left block
    ###########################
    cs = coupling(length,length+1) # list with couplings
    lS_plus_lB = (tensorial_operator(id_lS,lB) + 
                   tensorial_operator(lS,id_lB) + 
                   coupled_hamiltonian(lS_lB,lB_lS,cs=cs) ) # block+right site 
    lS_lB_new = [tensorial_operator(l,id_lB) for l in lS_lB] # right site to block 
    ###########################
    # now couple the two halves
    ###########################
# left and right parts
    cs = coupling(length+1,length+2) # list with couplings
    lSB_rSB = (tensorial_operator(lS_plus_lB,id_rSB,LO=LO) + 
                       tensorial_operator(id_lSB,rS_plus_rB,LO=LO) )
    # promote coupling between sites to the new dimension
    lS_rS_new = [tensorial_operator(r,id_lB) for r in lS_rS] # left site 
    rS_lS_new = [tensorial_operator(r,id_rB) for r in rS_lS] # left site 

    lS_plus_rS = coupled_hamiltonian(lS_rS_new,rS_lS_new,cs=cs,LO=LO) # left times right  
    lSB_rSB = lSB_rSB + lS_plus_rS # couple left and right sites
    outdict["hamiltonian"] = lSB_rSB # hamiltonian
    # ground state and project on DMRG manifold
    diml,dimr = dim_lB*dim_lS, dim_rB*dim_rS  # right and left dimensions
    gsout = ground_state(lSB_rSB,diml,dimr,v0=None,
                                   dmrgp=outdict) # get the smallest state
    eout = gsout.energy # GS energy
    # store several quantities
    dmrgout.iteration_energy.append(gsout.energy) # store energy
    dmrgout.energy = gsout.energy # store energy
    dmrgout.energies = gsout.energies # store energies
    dmrgout.iteration_energies.append(gsout.energies) # store energies
    dmrgout.iteration_dmdis.append(gsout.dmdis) # store DM distance
    #########################
    es = gsout.energies # excited states eenrgies
    wf0 = gsout.wf # wavefunction
    dmat = gsout.dm  # density matrix
    e0 = es[0] # ground state energy
    sisj = tensorial.exp_val(wf0,lS_plus_rS) # calculate the correlator
    if not silent:
      print("Correlation",sisj)
      print("Energy per site",e0/(dmrgp["total_length"]+2))
      print("Selected energy",eout)
      print()
    # now look for the relevant subspace
    R,entropy = states_dmat(dmat,dmrgp) # get the proj onto the main space
    dmrgout.iteration_entropies.append(entropy) # store energies
    # function to project matrices
    def old2new(A,mtype=None,iterative=False):
      """Project a certain operator into the new basis"""
      if iterative: return [old2new(Ai,mtype=mtype) for Ai in A]
      if mtype=="lSB": # left site+block
        out = A
      elif mtype=="lB": # left site+block
        out = tensorial_operator(id_lS,A)
      elif mtype=="lS": # left site+block
        out = tensorial_operator(A,id_lB)
      else: raise # not considered
      return project(out,R) # project the operator

    # put the right and left couplings in the full basis
    new_lB_lS = [old2new(M,mtype="lS") for M in lS_rS] # left block operators, projected
    new_lB = old2new(lS_plus_lB,mtype="lSB") # left block Hamiltonian, projected
    # store the results
    outdict["energy"] = eout
    outdict["energies"] = es - eout
#    outdict["v0"] = wf0 # save wavefunction
    outdict["correlator"] = sisj
    outdict["entropy"] = entropy
    outdict["lB_lS"] = new_lB_lS # coupling to the left block
    outdict["lB"] = new_lB # left hamiltonian
    outdict["length"] += 1 # the length if the chain has increased by 2 units
    # now store the site operators, in full form (lSB basis)
    if callable(dmrgp["target_function"]): # if there is a function
      opsSB = [old2new(oB,mtype="lB") + old2new(oS,mtype="lS") for (oB,oS) in 
                               zip(dmrgp["sum_block_operators"],dmrgp["left_site_operators"])]
      outdict["sum_block_operators"] = opsSB # store 
    # first transform the old ones
    if dmrgp["store_operators"]:
      outdict["density_matrix"] = old2new(dmat,mtype="lSB") # store dmat
      outdict["wf"] = wf0 # save wavefunction
      outdict["block_operators"] = [[old2new(m,mtype="lB") for m in s ] 
                                   for s in dmrgp["block_operators"]]
    # and add the new site
      ops = dmrgp["site_operator_generator"](0)
      outdict["site_operators"] = [old2new(op,mtype="lS") for op in ops]
      outdict["block_operators"].append(outdict["site_operators"])
#      print("# of stored operators",len(outdict["block_operators"]))
      if outdict["store_full_operators"]: # if full operators should be stored
        outdict["has_full_operators"] = True # has full operators
        outdict["left_block_operators_full"] = [[m2full(m,mtype="lB") for m in s ] 
                                   for s in dmrgp["block_operators"]]
        outdict["left_site_operators_full"] = [m2full(op,mtype="lS") for op in ops]
        outdict["store_full_operators"] = False # do not do it the next time
   
    return outdict # return the dictionary


#
#
#
#def dmrg_BoBo(dmrgp,integrate_right=True,LO=False):
#    """Perform a single DMRG step.
#    The geometry used in this function is
#            [ ] o [ ] o
#    using periodic boundary conditions"""
#    outdict = deepcopy(dmrgp) # output dictionary
#    ##### Couplings #####
#    rS_rB = dmrgp["rS_rB"] 
#    lS_lB = dmrgp["lS_lB"] 
#    lS_rB = dmrgp["lS_rB"] 
#    rB_rS = dmrgp["rB_rS"] 
#    rB_lS = dmrgp["rB_lS"] 
#    rS_lB = dmrgp["rS_lB"] 
#    lB_lS = dmrgp["lB_lS"] 
#    lB_rS = dmrgp["lB_rS"] 
#    #### Onsite Hamiltonians
#    rB = dmrgp["rB"] # right hamiltonian
#    lB = dmrgp["lB"] # left hamiltonian
#    rS = dmrgp["rS"] # left hamiltonian
#    lS = dmrgp["lS"] # left hamiltonian
#    nstates = dmrgp["number_of_states"] # number of states retained
#    target = dmrgp["target"] # number of target states
#    retain_states = dmrgp["retain_states"] # number of target states
#    # dimensions
#    dim_lS = lS.shape[0] # left site dimension
#    dim_rS = rS.shape[0] # right site dimension
#    dim_rB = rB.shape[0] # dimension of the block 
#    dim_lB = lB.shape[0] # dimension of the block 
#    # idensity operators
#    id_lB = sp.eye(dim_lB,dtype=np.complex) # identity operator
#    id_rB = sp.eye(dim_rB,dtype=np.complex) # identity operator
#    id_lS = sp.eye(dim_lS,dtype=np.complex) # identity operator
#    id_rS = sp.eye(dim_rS,dtype=np.complex) # identity operator
#    id_rSB = sp.eye(dim_rS*dim_rB,dtype=np.complex) # right identity operator
#    id_lSB = sp.eye(dim_lS*dim_lB,dtype=np.complex) # right identity operator
#    ###########################
#    # couple site to right block
#    ###########################
#    rS_plus_rB = (tensorial_operator(id_rS,rB) +
#                   tensorial_operator(rS,id_rB) + 
#                   coupled_hamiltonian(rS_rB,rB_rS) ) # block+right site 
#    ###########################
#    # couple site to left block
#    ###########################
#    lS_plus_lB = (tensorial_operator(id_lS,lB) + 
#                   tensorial_operator(lS,id_lB) + 
#                   coupled_hamiltonian(lS_lB,lB_lS) ) # block+right site 
#    ###########################
#    # now couple the two halves
#    ###########################
## left and right parts
#    LO = True
#    lSB_rSB = (tensorial_operator(lS_plus_lB,id_rSB,LO=LO) + 
#                       tensorial_operator(id_lSB,rS_plus_rB,LO=LO) )
#    # promote coupling between sites to the new dimension
#    lS_rB_new = [tensorial_operator(r,id_lB) for r in lS_rB]  
#    lB_rS_new = [tensorial_operator(id_lS,r) for r in lB_rS]  
#    rS_lB_new = [tensorial_operator(r,id_rB) for r in rS_lB]  
#    lB_rS_new = [tensorial_operator(id_lS,r) for r in lB_rS]  
#    rB_lS_new = [tensorial_operator(id_rS,r) for r in rB_lS]  
#
#    lS_plus_rB = coupled_hamiltonian(lS_rB_new,rB_lS_new,LO=LO)   
#    lB_plus_rS = coupled_hamiltonian(lB_rS_new,rS_lB_new,LO=LO)   
#    lSB_rSB = lSB_rSB + lS_plus_rB + lB_plus_rS # couple left and right parts
#    # ground state and project on DMRG manifold
#    diml,dimr = dim_lB*dim_lS, dim_rB*dim_rS  # right and left dimensions
#    gsout = ground_state(lSB_rSB,diml,dimr,v0=None) # get the smallest state
#    eout = gsout.energy
#    es = gsout.energies
#    wf0 = gsout.wf
#    dmat = gsout.dm
#    Ain = tensorial_operator(coupled_hamiltonian(lS_lB,lB_lS),id_rSB,LO=LO)
#    sisj = tensorial.exp_val(wf0,lS_plus_rB) # calculate the correlator
##    Ain = tensorial_operator(tensorial_operator(lS,id_lB),id_rSB)
##    sisj = spectrum.exp_val(wf0,Ain) # calculate the correlator
#    print("Correlation",sisj)
#    print("Energy per site",e0/(2*dmrgp["length"]+2))
#    # now look for the relevant subspace
#    R,entropy = states_dmat(dmat,nstates=nstates) # get the proj onto the main space
#    # function to project matrices
#    def old2new(A,mtype=None):
#      """Project a certain operator into the new basis"""
#      if mtype=="lSB": # left site+block
#        out = A
#      elif mtype=="lB": # left site+block
#        out = tensorial_operator(id_lS,A)
#      elif mtype=="lS": # left site+block
#        out = tensorial_operator(A,id_lB)
#      else: raise # not considered
#      return project(out,R) # project the operator
#
#    # put the couplings in the new basis
#    new_lB_rS = [old2new(M,mtype="lB") for M in lB_rS] # left block operators, projected
#    new_lB_lS = [old2new(M,mtype="lS") for M in lS_rB] # left block operators, projected
#    new_rB_lS = new_lB_rS
#    new_rB_rS = new_lB_lS 
#    new_lB = old2new(lS_plus_lB,mtype="lSB") # left block Hamiltonian, projected
#    new_rB = new_lB # right block Hamiltonian, projected
#    # store the results
#    outdict["energy"] = e0
#    outdict["correlator"] = sisj
#    outdict["entropy"] = entropy
#    outdict["lB_rS"] = new_lB_rS # coupling to the left block
#    outdict["lB_lS"] = new_lB_lS # coupling to the left block
#    outdict["rB_lS"] = new_rB_lS # coupling to the left block
#    outdict["rB_rS"] = new_rB_rS # coupling to the left block
#    outdict["rB"] = new_rB # right hamiltonian
#    outdict["lB"] = new_lB # left hamiltonian
#    outdict["length"] += 1 # the length if the chain has increased by 2 units
#    outdict["density_matrix"] = old2new(dmat,mtype="lSB") # store dmat in new basis
#    # now store the site operators, in full form (lSB basis)
#    # first transform the old ones
#    outdict["site_operators"] = [[old2new(m,mtype="lB") for m in s ] 
#                                   for s in dmrgp["site_operators"]]
#    # and add the new site
#    ops = dmrgp["site_operator_generator"](0)
#    outdict["site_operators"].append([old2new(op,mtype="lS") for op in ops])  
#    print("Stored sites",len(outdict["site_operators"]))
#    return outdict # return the dictionary
#
#
#





def initial_BooB(right,left,onsite,
        site_operators=lambda i: [],cs=None):
  """Generate the initial dictionary"""
  if cs is None: raise
  rS_rB = right(0)
  rS_lS = left(0)
  lS_rS = right(0)
  lS_lB = left(0)
  # identity operators
  dim_rS = rS_rB[0].shape[0] # dimension of right subspace
  dim_lS = lS_rS[0].shape[0] # dimension of left subspace
  id_rS = sp.eye(dim_rS,dtype=np.complex) # left identity operator
  id_lS = sp.eye(dim_lS,dtype=np.complex) # left identity operator
  # initial Hamiltonian
  Block = (coupled_hamiltonian(right(0),left(1),cs=cs) +
            tensorial_operator(id_lS,onsite(1)) +
            tensorial_operator(onsite(0),id_rS) )
  #################################################
  # create the first iteration of block couplings #
  #################################################
  lB_lS = [tensorial_operator(id_lS,r) for r in right(1)]
  rB_rS = lB_lS # same
  dmrgp = dict()
  dmrgp["lS_lB"] = left(2) # coupling to the left block
  dmrgp["lB_lS"] = lB_lS # coupling to the left block
  dmrgp["lB"] = Block.copy() # left block hamiltonian
  dmrgp["lS"] = onsite(2) # left site hamiltonian
  dmrgp["v0"] = None # initial wavefunction
# site operators 
  ops0 = [tensorial_operator(O,id_rS) for O in site_operators(0)]
  ops1 = [tensorial_operator(id_lS,O) for O in site_operators(1)]
  dmrgp["block_operators"] = [ops0,ops1] 
  dmrgp["sum_block_operators"] = [o0+o1 for (o0,o1) in zip(ops0,ops1)]
  # keep in mind that the right site is actually the number 0,
  # whereas the left one is number 1, because the initial geometry is
  # assumed o [01], and this will grow form the left
  # generating function
#  dmrgp["site_operator_generator"] = site_operators 
  dmrgp["store_operators"] = True 
  dmrgp["length"] = 2
  dmrgp["store_full_operators"] = False # store the full operators
  dmrgp["has_full_operators"] = False # has the full operators
  return dmrgp




def initial_BoBo(right,left,onsite,site_operators=lambda i: []):
  """Generate the initial dictionary"""
  rS_rB = right(0)
  rS_lS = left(0)
  lS_rS = right(0)
  lS_lB = left(0)
  # identity operators
  dim_rS = rS_rB[0].shape[0] # dimension of right subspace
  dim_lS = lS_rS[0].shape[0] # dimension of left subspace
  id_rS = sp.eye(dim_rS,dtype=np.complex) # left identity operator
  id_lS = sp.eye(dim_lS,dtype=np.complex) # left identity operator
  # insitial Hamiltonian
  Block = (coupled_hamiltonian(lS_rS,rS_lS) +
            tensorial_operator(id_lS,onsite(0)) +
            tensorial_operator(onsite(0),id_rS) )
  #################################################
  # create the first iteration of block couplings #
  #################################################
  rB_rS = [tensorial_operator(l,id_rS) for l in rS_lS]
  lB_lS = rB_rS # same
  rB_lS = [tensorial_operator(id_lS,r) for r in lS_rS]
  lB_rS = rB_lS # same
  lS_rB = left(0) # coupling to the right
  rS_lB = lS_rB # coupling to the right
  wf0 = None # no initial vector
  rB = Block.copy()
  lB = Block.copy()
  rS = onsite(0)
  lS = onsite(0)
  dmrgp = dict()
  dmrgp["rS_rB"] = rS_rB # coupling to the right bloxk
  dmrgp["lS_rB"] = lS_rB # coupling between sites
  dmrgp["rS_lB"] = rS_lB # coupling between sites
  dmrgp["lS_lB"] = lS_lB # coupling to the left block
  dmrgp["rB_rS"] = rB_rS # coupling to the right bloxk
  dmrgp["rB_lS"] = rB_lS # coupling to the right bloxk
  dmrgp["lB_lS"] = lB_lS # coupling to the left block
  dmrgp["lB_rS"] = lB_rS # coupling to the left block
  dmrgp["rB"] = rB # right block hamiltonian
  dmrgp["lB"] = lB # left block hamiltonian
  dmrgp["rS"] = rS # right site hamiltonian
  dmrgp["lS"] = lS # left site hamiltonian
  # site operators
  opsr = [tensorial_operator(id_lS,O) for O in site_operators(1)]
  opsl = [tensorial_operator(O,id_rS) for O in site_operators(0)]
  dmrgp["site_operators"] = [opsr,opsl] 
  dmrgp["site_operator_generator"] = site_operators 
  dmrgp["store_operators"] = False # do not store the operators
  dmrgp["store_full_operators"] = False # store the full operators
  dmrgp["has_full_operators"] = False # has the full operators
  return dmrgp




