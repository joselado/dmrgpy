import numpy as np
from . import operatornames


def get_cached_excited_states(self,n=20,scale=10.0,**kwargs):
    """Wrapper around self.get_excited_states(purify=False,...) that
    reuses a previous excited-state search on the same chain instance
    when the (n,scale,gram_schmidt) settings match. The excited states
    themselves don't depend on the two operators A,B being correlated,
    so computing them again for every dynamical_correlator(name=(A,B))
    call on the same Hamiltonian -- e.g. looping over several operator
    pairs, as several examples/dynamical_correlator/*_excited scripts do
    -- used to rerun the whole O(n) sequential DMRG excited-state search
    from scratch each time. The cache is invalidated by
    Many_Body_Chain.restart()/set_hamiltonian (see manybodychain.py),
    since a new ground state changes what "the excited states" means."""
    key = (n,scale,getattr(self,"excited_gram_schmidt",False))
    cache = getattr(self,"_dcex_excited_cache",None)
    if cache is not None and cache[0]==key:
        return cache[1]
    result = self.get_excited_states(n=n,purify=False,scale=scale,**kwargs)
    self._dcex_excited_cache = (key,result)
    return result


def dynamical_correlator(self,name="XX",i=0,j=0,delta=2e-2,
        nex=20,es=np.linspace(-1.0,10.0,1000),scale=10.0,**kwargs):
    """
    Compute dynamical correlator with excited states with DMRG
    """
#    self.gs_energy()
#    name = operatornames.str2MO(self,name)
#    task = {"dynamical_correlator_excited":"true",
#            "scale_lagrange_excited":str(scale),
#            "nexcited":str(nex),
#            }
#    self.task = task # override tasks
#    name[0] = name[0].get_dagger()
#    self.execute(lambda: name[0].write(name="dc_multioperator_i.in"))
#    self.execute(lambda: name[1].write(name="dc_multioperator_j.in"))
#    self.execute( lambda : taskdmrg.write_tasks(self)) # write tasks
#    self.execute( lambda : self.run()) # run calculation
#    # now read the data
#    eex = self.get_file("EXCITED.OUT").T[0] # excitation energies
#    eex = eex[1:len(eex)] - eex[0] # substract GS energy
#    cs = self.get_file("EXCITED_OVERLAPS.OUT") # overlaps with GS
#    c1 = cs[:,0] + 1j*cs[:,1] # first overlap
#    c2 = cs[:,2] - 1j*cs[:,3] # second overlap
    # compute the two correlators. "scale" here used to be silently
    # dropped (never forwarded to get_excited_states, due to this
    # function's own "scale" kwarg shadowing it) -- now passed through
    # explicitly, with the default kept at 10.0 to match what
    # get_excited_states_dmrg's own default actually was, so this fix
    # doesn't change default behavior, only lets a caller-supplied
    # scale= finally take effect.
    esex,wsex = get_cached_excited_states(self,n=nex,scale=scale,**kwargs)
    # normalize name= first: a documented string ("ZZ" with i=/j=) used to
    # reach name[0].get_dagger() as a bare str and die with
    # "AttributeError: 'str' object has no attribute 'get_dagger'"
    name = operatornames.str2MO(self,name,i=i,j=j)
    A,B = name[0].get_dagger(),name[1] # operators
    h = self.hamiltonian # Hamiltonian
    n = len(wsex) # number of states
    # Matrix representations in the raw nex-state basis returned by the
    # excited-state search, WITHOUT assuming it is orthonormal: DMRG
    # excited states are only approximately orthogonal even after
    # Gram-Schmidt (each MPS is separately bond-truncated). wfi.aMb(Op,wfj)
    # (a direct <i|Op|j> sandwich contraction) replaces what used to be
    # wfi.dot(Op*wfj): the latter first applies Op to wfj as an
    # intermediate MPS, requiring a variational compression/fit at bond
    # dimension maxm for every one of the nex^2 matrix elements -- aMb
    # needs no such intermediate MPS (and thus no compression error) for
    # any of Sop/Aop/Bop/Hop.
    Sop = np.array([[wsex[a].dot(wsex[b]) for b in range(n)] for a in range(n)])
    Aop = np.array([[wsex[a].aMb(A,wsex[b]) for b in range(n)] for a in range(n)])
    Bop = np.array([[wsex[a].aMb(B,wsex[b]) for b in range(n)] for a in range(n)])
    Hop = np.array([[wsex[a].aMb(h,wsex[b]) for b in range(n)] for a in range(n)])
    # H is diagonalized in this (possibly non-orthogonal) basis via the
    # generalized eigenvalue problem Hc=eSc, solved by canonical (Loewdin)
    # orthogonalization rather than scipy.linalg.eigh(Hop,b=Sop) directly:
    # that call needs a Cholesky factorization of Sop and either raises
    # LinAlgError ("leading minor ... not positive definite") or silently
    # returns a spurious huge-magnitude eigenvalue when Sop is
    # (near-)singular -- confirmed directly, e.g. a 2-state basis with
    # overlap 0.999999999 already produces a spurious eigenvalue of
    # ~1.3e9. Near-linear-dependence among the nex excited states is
    # expected in practice (nex requested past what the bond dimension can
    # resolve as distinct, or a genuine physical degeneracy), so this isn't
    # a hypothetical edge case: it's exactly the situation the old
    # gram_smith()+remove_none() combination used to guard against, by
    # dropping near-zero-norm directions during orthogonalization.
    # Diagonalizing Sop and keeping only directions above a relative
    # tolerance reproduces that same protection: U's columns are an
    # explicitly S-orthonormal (U^H Sop U = 1) basis for the
    # well-conditioned part of the original nex-dim space, so Hop
    # transforms into a smaller, genuinely orthonormal basis where a plain
    # (non-generalized) eigh is exact and numerically safe.
    from scipy.linalg import eigh
    sval,svec = eigh(Sop)
    keep = sval>1e-8*np.max(sval)
    U = svec[:,keep]/np.sqrt(sval[keep])
    Hred = np.conjugate(U.T)@Hop@U
    esex,vs = eigh(Hred) # well-conditioned, ordinary eigenvalue problem
    wsex = (U@vs).T # coefficients back in the original nex-state basis
    # from now on we operate with numpy arrays
    wf0 = wsex[0]
    wfa = Aop@wf0 # A times ground state
    wfb = Bop@wf0 # B times ground state
    c1 = [np.conjugate(wfa).dot(wfi) for wfi in wsex] # matrix element
    c2 = [np.conjugate(wfi).dot(wfb) for wfi in wsex] # matrix element
    eex = esex - esex[0] # difference
    es,adv = dcex(eex,c1,c2,es=es,delta=delta) # return correlator
    es,ret = dcex(eex,c2,c1,es=es,delta=-delta) # return correlator
    return es,1j*(adv-ret)/(2.*np.pi) # return correlator





def dcex(eex,c1,c2,es=np.linspace(-1.0,10.0,300),delta=1e-1):
    """
    Compute a dynamical correlator with excitation energies and
    matrix elements
    """
    if len(eex)!=len(c1): raise # something wrong
    if len(eex)!=len(c2): raise # something wrong
    if es is None: es = np.linspace(-1.0,np.max(eex),2000)
    out = np.zeros(len(es),dtype=np.complex128) # output
    out = [] # empty list
    out = 0.0j # initialize
    for i in range(len(eex)): # loop over excited states
        out = out + c1[i]*c2[i]/(es-eex[i]+1j*delta) # dynamical correlator
    return es,out # return 








