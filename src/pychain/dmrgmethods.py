
def dmrg_BooB(indict,integrate_right=True):
    """Perform a single DMRG step"""
    outdict = deepcopy(indict) # output dictionary
    rs = indict["right_site_coupling"] # coupling in the right site
    ls = indict["left_site_coupling"] # coupling in the left site
    crs = indict["right_block_coupling"] # coupling to the right bloxk
    cls = indict["left_block_coupling"] # coupling to the left block
    hr = indict["right_block"] # right hamiltonian
    hl = indict["left_block"] # left hamiltonian
    rons = indict["right_site"] # left hamiltonian
    lons = indict["left_site"] # left hamiltonian
    nstates = indict["number_of_states"] # number of states retained
    target = indict["target_states"] # number of target states
    # perform the calculation
    diml = ls[0].shape[0] # left site dimension
    dimr = rs[0].shape[0] # right site dimension
    dimc = hr.shape[0] # dimension of the block 
    idc = sp.eye(hr.shape[0],dtype=np.complex) # block identity operator
    idl = sp.eye(diml,dtype=np.complex) # block identity operator
    idr = sp.eye(dimr,dtype=np.complex) # block identity operator
    # couple site to right block
    hr = tensorial_operator(idr,hr) # block+right site 
    if not rons is None: # add right onsite
      hr = hr + tensorial_operator(rons,idc) # onsite of the new rigth site
    hr = hr + coupled_hamiltonian(rs,crs) # couple right site to block
    rs = [tensorial_operator(r,idc) for r in rs] # right site to block 
    # couple site to left block
    hl = tensorial_operator(idl,hl) # block+left site 
    if not lons is None: # add left onsite
      hl = hl + tensorial_operator(lons,idc) # onsite of the new left site
    hl = hl + coupled_hamiltonian(ls,cls) # couple site to block
    ls = [tensorial_operator(l,idc) for l in ls] # left site to block 
    # now couple the two halves
    idr2 = sp.eye(hr.shape[0],dtype=np.complex) # right identity operator
    idl2 = sp.eye(hl.shape[0],dtype=np.complex) # right identity operator
# left and right parts
    h = tensorial_operator(hl,idr2) + tensorial_operator(idl2,hr)
#    print(np.sum(hl-hr))
    htwo = coupled_hamiltonian(ls,rs) # left times right  
    h = h + htwo # couple left and right sites
    # ground state and project on DMRG manifold
    e0,wf0 = ground_state(h,v0=None,target=target) # get the smallest state
    sisj = spectrum.exp_val(wf0,csc_matrix(htwo)) # calculate the correlator
    dmat = trace_over(wf0,dimc*diml,dimc*dimr) # trace over right site + block
    entropy = get_entropy(dmat)
    # now look for the relevant subspace
    R = states_dmat(dmat,nstates=nstates) # get the proj onto the main space
    # put the right and left couplings in the full basis
    new_crs = [project(c,R) for c in rs] # right block operators, projected
    new_hr = project(hr,R) # right block Hamiltonian, projected
    new_hl = project(hl,R) # left block Hamiltonian, projected
    new_cls = [project(c,R) for c in ls] # left block operators, projected
    # store the results
    outdict["energy"] = e0
    outdict["correlator"] = sisj
    outdict["entropy"] = entropy
    outdict["right_block_coupling"] = new_crs # coupling to the right bloxk
    outdict["left_block_coupling"] = new_cls # coupling to the left block
    outdict["right_block"] = new_hr # right hamiltonian
    outdict["left_block"] = new_hl # left hamiltonian
    outdict["right_site"] = rons # right site hamiltonian
    outdict["left_site"] = lons # left site hamiltonian
    outdict["length"] += 2 # the length if the chain has increased by 2 units
    # function to express a certain operator in the new basis
    def old2new(A,site="left_block"):
      """Project a certain operator into the new basis"""
      if site=="left_block": # for right blocks
        B = tensorial_operator(idl,A)
        Id = tensorial_operator(idr,idc)
        B = tensorial_operator(B,Id) # operator in new basis
      else: raise # not considered
      return project(B,R) # project the operator
    return outdict # return the dictionary





