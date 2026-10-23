import numpy as np



# dynamical_correlator_kpm() used to live here. It read a KPM_ENTROPY.OUT
# file written by the old file-based backend, via self.execute() -- both of
# which were removed with that backend, so it could only ever raise
# AttributeError. Its one caller (examples/dynamical_correlator/
# dynamical_correlator_entropy) no longer uses it.

# compute_entropy(b=None), a geometric mean over the bonds, lived here. It
# had no caller, and disagreed with the arithmetic mean MPS.get_entropy()
# returns (0 as soon as one bond was unentangled), so it was removed.


def compute_entropy_single(self,psi,b=1):
    """Compute entanglement entropy in a bond"""
    # Bond b is the cut between sites b-1 and b, so with sites indexed
    # 0..ns-1 the valid range is 1..ns-1 -- which is exactly what
    # MPS.get_entropy(b=None) loops over (bonds 1..ns-1). The guard
    # used to read b>self.ns, one too far: b==ns slipped through to
    # ITensor, whose own check calls abort(), killing the whole process
    # with an uncatchable SIGABRT instead of reporting a bad index (the
    # "python" backend raised a plain IndexError for the same call).
    if b<1 or b>=self.ns:
        raise IndexError("bond %s is out of range: this chain has %d "
                "sites, so its bonds are 1..%d (bond b is the cut "
                "between sites b-1 and b)"%(repr(b),self.ns,self.ns-1))
    from .mode import resolve_mode
    if resolve_mode(self)=="ED":
        # MPS-only, like the rest of this function: there is no ED
        # implementation of the bond entropy (the ED object has no MPS to
        # cut). Say so, instead of handing an ED State to
        # self._session and failing with the opaque "'State' object has
        # no attribute 'cpp_handle'" -- same treatment as
        # Many_Body_Chain.get_distribution_moments. Ahead of the
        # julia_live branch, which used to come first and so failed the
        # same way on 'jlmps' (2026-09-25 hole hunt 25b, lead
        # session-julia-ed-guard-carveout-rdm-bond-entropy).
        raise NotImplementedError(
            "the bond entanglement entropy has no ED implementation (it "
            "cuts an MPS bond, and the ED backend has no MPS). Note "
            "mode.py routes to ED on its own when the requested C++ "
            "extension is unavailable, or for itensor_version=3 on a "
            "chain with fewer than 3 sites. Use get_site_entropy/"
            "get_pair_entropy, which do have an ED route.")
    if self.itensor_version=="julia_live":
        from .mpsjulialive import entropy as entjl
        return np.abs(entjl.bond_entropy(psi,b))
    # The entropy belongs to the ray, so cut the normalized state: v2 and v3
    # take -sum w ln w over the raw Schmidt weights, and a state of norm^2 4
    # came out up to 3.5 too high there while "python" (which normalizes)
    # was exact (2026-10-23 hole hunt, finding 4). A unit-norm state, every
    # solver output, is left untouched.
    n2 = psi.dot(psi).real
    if not n2>0.: raise ValueError("bond entropy of a zero state")
    if abs(n2-1.)>1e-12: psi = psi*(1./np.sqrt(n2))
    return np.abs(self._session.bond_entropy(psi.cpp_handle,b))



def bond_entropy(self,wf,i,j):
    """Compute the entropy of a state in bond i"""
    # validate the site indices the caller actually passed, before they
    # are collapsed into a single bond index by max() -- otherwise the
    # error message points at a bond the caller never named
    for k in (i,j):
        if not (0<=k<self.ns):
            raise IndexError("site %s is out of range: this chain has %d "
                    "sites, indexed 0..%d"%(repr(k),self.ns,self.ns-1))
    if abs(i-j)==1: # use the DMRG approach
        return compute_entropy_single(self,wf,b=max([i,j]))
    else: raise ValueError("get_bond_entropy needs two adjacent sites "
            "(|i-j|==1), got i=%s, j=%s"%(repr(i),repr(j)))
#    from .densitymatrix import reduced_dm_projective
#    dm = reduced_dm_projective(self,wf,i=i,j=j) # compute density matrix
#    return entropy_dm(dm,normalize=True) # return the entropy


def _unit_state(wf):
    """The state divided by its norm; the entropy belongs to the ray"""
    n2 = wf.dot(wf).real
    if not n2>0.: raise ValueError("entropy of a zero state")
    if abs(n2-1.)>1e-12: wf = wf*(1./np.sqrt(n2))
    return wf


def pair_entropy(self,wf,i,j):
    """Compute the entropy of a state in bond i"""
    from .densitymatrix import reduced_dm_projective
    # normalized like site_entropy: a state of norm^2 4 raised, and one of
    # norm^2 1.0004 passed the trace check and came out 1.1e-4 too high
    wf = _unit_state(wf)
    dm = reduced_dm_projective(self,wf,i=i,j=j) # compute density matrix
    return entropy_dm(dm,normalize=False) # return the entropy


def site_entropy(self,wf,i):
    """Compute the entropy of a state in bond i"""
    from .densitymatrix import reduced_dm_projective
    wf = _unit_state(wf) # normalize the wavefunction
    dm = reduced_dm_projective(self,wf,i=i,j=None) # compute density matrix
    return entropy_dm(dm,normalize=False) # return the entropy


def mutual_information(self,wf,i,j):
    """Compute the mutual information"""
    si = site_entropy(self,wf,i) # entropy in i
    sj = site_entropy(self,wf,j) # entropy in j
    sij = pair_entropy(self,wf,i,j) # joint entropy
    return si + sj - sij # return the mutual information



def entropy_dm(dm,normalize=False):
    from scipy.linalg import eigvalsh
    if np.abs(1.-np.trace(dm))>1e-3:
        raise ValueError("this density matrix is not normalized (trace = "
                +str(np.trace(dm))+", expected 1); its entropy is not "
                "meaningful")
    ds = eigvalsh(dm) # compute eigenvalues
    # Relative floor, only to drop roundoff (including negative eigenvalues).
    # An absolute 1e-6 dropped the -p ln p of a real weight p<=1e-6, so a
    # nearly product site with p=5.6e-7 gave 5.6e-7 for an exact 8.6e-6.
    ds = ds[ds>1e-14*np.max(np.abs(ds))]
    if normalize: 
        n = dm.shape[0]
        norm = np.log(n) # normalize
    else: norm = 1.
    return -np.sum(ds*np.log(ds))/norm # return the entropy


def central_charge(wf):
    """Compute the central charge of a wavefunction assuming it is critical"""
    L = len(wf.MBO.sites) # number of sites
    sr = [wf.get_bond_entropy(i-1,i) for i in range(1,L)] # entropy
    sr = np.array(sr)
    ls = np.array(range(1,L)) # lengths
    # S = c/6 log(2L/pi sin(pi l/L)) + const, J.Stat.Mech.0406:P06002,2004.
    # Linear in (c, const), so solved as a linear least-squares problem,
    # deterministic and with the same minimum the old random-start minimizer
    # found (identical to 6 digits on XX chains of 8, 20 and 40 sites). On a
    # short chain the alternating even/odd term dominates and that minimum
    # can be negative: -0.63 on the 8-site Heisenberg chain.
    x = np.log(2*L/np.pi*np.sin(np.pi*ls/L))/6.
    A = np.array([x,np.ones(len(x))]).T
    return np.linalg.lstsq(A,sr,rcond=None)[0][0] # return central charge

