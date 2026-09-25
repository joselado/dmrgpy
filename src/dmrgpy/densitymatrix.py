# routines to compute density matrices
import numpy as np

def reduced_dm(self,i=0,mode="DMRG"):
    """
    Compute the reduced density matrix
    """
    from .mode import resolve_mode
    # On every backend, julia_live included: it used to be exempted
    # (`self.itensor_version!="julia_live" and ...`), so on julia_live
    # sc.mode="ED" died with "'State' object has no attribute 'jlmps'"
    # and get_rdm(mode="ED") quietly returned the DMRG matrix (2026-09-25
    # hole hunt 25b, lead session-julia-ed-guard-carveout-rdm-bond-
    # entropy). resolve_mode() answers "DMRG" on a julia_live chain that
    # asks for nothing else, so this costs that backend nothing.
    if resolve_mode(self,mode=mode)=="ED":
        # No ED implementation exists (grep edtk/: there is no
        # reduced_dm there), and this function is session-only, so
        # without this guard an ED State reached self._session and died
        # with the opaque "'State' object has no attribute 'cpp_handle'"
        # -- the very symptom get_distribution's own get_mode() fix
        # quotes. Say what happened instead, in that method's wording.
        # Note reduced_dm_projective() below *does* work under ED and is
        # what get_site_entropy/get_pair_entropy use, but it is not a
        # drop-in: its basis is the projector list ([N,Cdag] /
        # [Sz+1/2,S+]), so its diagonal comes out in the opposite order
        # from ITensor's (empty,occupied) convention.
        raise NotImplementedError(
            "get_rdm has no ED implementation (the reduced density "
            "matrix is read off an MPS bond, and the ED backend has no "
            "MPS). Note mode.py routes to ED on its own when the "
            "requested C++ extension is unavailable, or for "
            "itensor_version=3 on a chain with fewer than 3 sites.")
    wf = self.get_gs() # compute ground state
    # A density matrix of the ray: every backend's reduced_dm divided the
    # state by <wf|wf> rather than by its square root, a no-op only while
    # the state is unit norm, and set_gs(c*s)/set_initial_wf(c*s)/
    # gs_energy(wf0=c*s, reconverge=False) hand over a caller's state
    # unswept, so the matrix came back as rho/c^2 (2026-09-25 hole hunt
    # 25b, finding 11). Normalized here once, for all four backends; the
    # backend's own division then divides by 1. Not wf.normalize(), which
    # returns None below its tolerance.
    nrm2 = float(np.real(wf.dot(wf)))
    if not nrm2>0.0: # "not >", so that a NaN is caught too
        raise ValueError("get_rdm: the chain's state has zero norm "
                "(<wf|wf> = %r), so it has no density matrix"%nrm2)
    if abs(nrm2-1.0)>1e-14: wf = wf*(1.0/np.sqrt(nrm2))
    if self.itensor_version=="julia_live":
        from .mpsjulialive import densitymatrix as dmjl
        return dmjl.reduced_dm(self,wf,i)
    return self._session.reduced_dm(wf.cpp_handle,i+1)



def reduced_dm_projective(self,wf,i=0,j=None):
    """Compute the reduced density matrix using a brute force approach"""
    from .fermionchain import Fermionic_Chain
    from .fermionchain import Spinful_Fermionic_Chain
    from .spinchain import Spin_Chain
    def projectors(k):
        """Build the projectors onto the different single-site components"""
        site = self.sites[k] # take this site
        if type(self)==Spin_Chain: # spin chain object
          if site==2: # S=1/2 site
              Szk = self.Sz[k] # Sz operator
              P01 = self.Sx[k]+1j*self.Sy[k] # project and rotate
              # we need to project on up/dn and rotate to up!
              return [Szk+0.5,P01] 
          elif site==3: # S=1 site
              raise NotImplementedError("reduced_dm_projective: S=1 "
                      "site projectors are not finished")
              Szk = self.Sz[k]
              return [Szk*(Szk+1)/2.,-(Szk-1)*(Szk+1),Szk*(1-Szk)/2.]
          else: raise NotImplementedError("reduced_dm_projective: no "
                  "projectors for site type "+str(site))
        elif type(self)==Fermionic_Chain: # spin chain object
          N = self.N[k] # density operator
          P01 = self.Cdag[k] # create an electron
          # we need to project on up/dn and rotate to up!
          return [N,P01] # return the projectors
        elif type(self)==Spinful_Fermionic_Chain: # spin chain object
          raise NotImplementedError("reduced_dm_projective: not "
                  "implemented yet for Spinful_Fermionic_Chain")
        else: raise NotImplementedError("reduced_dm_projective: not "
                "implemented for chain type "+type(self).__name__)
    Pi = projectors(i) # projectors for site i
    if j is not None: Pj = projectors(j) # projectors for site j
    else: Pj = [1] # workaround for a single site
    Pk = [] # projectors in the ij subspace
    for pi in Pi:
        for pj in Pj: Pk.append(pi*pj) # store the projector in this subspace
    # now compute the density matrix in this subspace
    n = len(Pk) # number of components
    dm = np.zeros((n,n),dtype=np.complex128) # initialize
    # dm[a,b] = <wf|Pa^dag Pb|wf>, evaluated as a single sandwich rather
    # than by first building the intermediate MPS Pa|wf>. Two reasons, one
    # of them load-bearing: the intermediate needs a variational
    # compression per element (avoidable error, the same argument
    # dcex.py's own aMb rewrite makes), and -- the actual bug -- some of
    # these projectors change the conserved charge (P01 = Sx+i*Sy is S+,
    # Cdag adds a particle), so applying one to the state made every
    # site/pair entropy and mutual information raise outright in
    # conserved-sector mode.
    in_sector = bool(getattr(self,"conserved_sector",None))
    for a in range(n):
        for b in range(n):
            op = Pk[a].get_dagger()*Pk[b]
            if in_sector:
                try: dm[a,b] = wf.aMb(op,wf)
                except (ValueError,NotImplementedError):
                    # Pa^dag Pb changes the sector's charge, so the chain
                    # cannot represent it -- and its expectation value in
                    # a fixed-charge state is exactly zero anyway, which
                    # is what that entry is
                    dm[a,b] = 0.0
            else:
                dm[a,b] = wf.aMb(op,wf)
#    print(dm.real)
    return dm # return density matrix



#
#def explicit_dm(sc,wf,inds=[0]):
#    """Compute the density matrix explicitly by summing
#    over all the vectors. This is a very heavy procedure, but good
#    for benchmarking and debugging"""
#    sc = sc.copy() # make a copy
#    for site in sc.sites: # check that you only have S=1/2
#        if site !=2: 
#            print("Only implemented for S=1/2")
#            raise # stop
#    # now loop over all the sites
#    for ii in inds: # loop over sites where you want the entropy
#        for Bz in [-1,1]: # the two combinations
#            Hi = sc.Sz[ii] # the two magnetic fields
#
#
#

