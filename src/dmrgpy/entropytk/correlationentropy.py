import numpy as np
import scipy.linalg as lg



def get_correlation_matrix(self,T=0.,**kwargs):
    """Compute the correlation matrix of a finite temperature state"""
    if T==0.: return get_correlation_matrix_zeroT(self,**kwargs)
    else: return get_correlation_matrix_finiteT(self,T=T,**kwargs)



def get_correlation_matrix_finiteT(self,T=1.,n=None,**kwargs):
    """Correlation matrix at finite temperature, ED only: sum
    Boltzmann-weighted per-eigenstate correlation matrices
    (get_correlation_matrix_zeroT with wf=eigenstate) over a set of
    exact eigenstates, exactly the same excited-states pattern as
    vevtk.thermalvev.thermal_vev_ex -- `n` controls how many eigenstates
    are summed over (default: the full Hilbert space whenever that's
    small enough to diagonalize exactly, matching algebra.maxsize;
    otherwise the same truncated n=10 that algebra.lowest_states falls
    back to on its own). `**kwargs` (operators=, basis=, dmmode=, ...)
    is forwarded to get_correlation_matrix_zeroT for each eigenstate,
    not to the excited-state search.

    A truncated sum whose highest-included state still carries
    non-negligible Boltzmann weight raises RuntimeError instead of
    silently returning a wrong average -- pass a larger n= in that case."""
    from ..algebra.algebra import maxsize
    dim = self.get_ED_obj().get_hamiltonian().shape[0] # full Hilbert space dimension
    if n is None: # caller did not request a specific truncation
        n_used = dim if dim<=maxsize else 10 # matches lowest_states' own default
    else: n_used = n
    truncated = n_used<dim # did we actually leave out any eigenstates?
    (es,wfs) = self.get_excited_states(mode="ED",n=n_used)
    es = es - np.min(es) # with respect to the GS
    beta = 1./T
    P = np.exp(-beta*es) # Boltzmann probabilities
    P = P/np.sum(P) # normalize
    if truncated and P[-1]>1e-6*np.max(P):
        raise RuntimeError(
            "get_correlation_matrix_finiteT: the "+str(n_used)+" excited states "
            "computed do not capture the full Boltzmann weight at T="+str(T)+
            " (highest included state has relative weight "+str(P[-1]/np.max(P))+
            " of the largest). Pass a larger n= to include more excited states.")
    dms = [get_correlation_matrix_zeroT(self,wf=wf,**kwargs) for wf in wfs]
    dm = np.zeros_like(dms[0])
    for i in range(len(es)): dm = dm + P[i]*dms[i] # Boltzmann-weighted average
    return dm # return matrix




def get_correlation_matrix_zeroT(self,operators=None,
                              basis ="electron", 
                              dmmode="fast",
                              wf=None,**kwargs):
    """Compute the correlation matrix of a ground state"""
    from .. import fermionchain
    if dmmode=="full" and basis=="Nambu":
        dmmode="fast"
        print("C++ mode not implemented with Nambu basis")
#    print(dmmode)
    if wf is None: wf = self.get_gs(**kwargs) # compute ground state
    wf = wf.normalize() # normalize wavefunction
    if operators is None: # no operators provided
        if fermionchain.isfermion(self):
            if basis=="Nambu": 
              operators = [o for o in self.C] 
              operators += [o for o in self.Cdag] 
            else: # just normal basis
              operators = self.C # fermionic operators
        else: 
            print("Unrecognized type",type(self))
            raise
    # create the matrix
    if dmmode=="simple":
        return correlation_matrix_clean(operators,wf,self)
    elif dmmode=="fast":
        return correlation_matrix_fast(operators,wf)
    elif dmmode=="explicit":
        return correlation_matrix_explicit(operators,wf)
    elif dmmode=="full":
        return cpp_correlation_matrix(wf)
    else: 
        print(dmmode,"not recognized")
        raise # not implemented
    n = len(operators)
    cm = np.zeros((n,n),dtype=np.complex128)
    for i in range(n):
        A = self.get_dagger(operators[i])
        for j in range(i,n):
            B = operators[j]
            wf1 = B*wf # first operator
            wf2 = A*wf1 # second operator
            out = wf.dot(wf2) # overlap
            cm[i,j] = out
            cm[j,i] = np.conjugate(out)
    return cm # return matrix



def correlation_matrix_explicit(operators,wf):
    """Compute the correlation entropy using a version that works
    for Julia ITensor"""
    n = len(operators) # number of operators
    As = operators # operators
    Bs = [A.get_dagger() for A in As] # dagger operators
    m = np.zeros((n,n),dtype=np.complex128) # initialize
    for i in range(n):
        for j in range(i,n):
            C = Bs[i]*As[j] # operator
            cij = wf.aMb(C,wf) # overlap
            m[i,j] = cij # original
            m[j,i] = np.conjugate(cij) # conjugate
    return m



def correlation_matrix_clean(operators,wf,self):
    """Compute the correlation matrix of a wavefunction with the fastest 
    algorithm"""
    # create the matrix
    n = len(operators)
    cm = np.zeros((n,n),dtype=np.complex128)
    for i in range(n):
        A = self.get_dagger(operators[i])
        for j in range(i,n):
            B = operators[j]
            wf1 = B*wf # first operator
            wf2 = A*wf1 # second operator
            out = wf.dot(wf2) # overlap
            cm[i,j] = out
            cm[j,i] = np.conjugate(out)
    return cm # return matrix

def correlation_matrix_fast(operators,wf):
    """Compute the correlation matrix of a wavefunction with the fastest 
    algorithm"""
    # create the matrix
    n = len(operators)
    wfs = [o*wf for o in operators] # compute all the wavefunctions
    cm = np.zeros((n,n),dtype=np.complex128)
    for i in range(n):
        for j in range(i,n):
            out = wfs[i].dot(wfs[j]) # overlap
            cm[i,j] = out
            cm[j,i] = np.conjugate(out)
    return cm # return matrix



## High order correlation entropy


def get_highorder_correlation_matrix(self,operators=None,wf=None,**kwargs):
    """Compute the correlation matrix of a ground state"""
    from .. import fermionchain
    if wf is None: wf = self.get_gs(**kwargs) # compute ground state
    if operators is None:
        if type(self)==fermionchain.Fermionic_Chain:
            operators = self.C # fermionic operators
        elif type(self)==fermionchain.Spinful_Fermionic_Chain:
            operators = self.C # fermionic operators
        else: raise
    # create the matrix
    n = len(operators)
    cm = np.zeros((n,n,n,n),dtype=np.complex128)
    for i in range(n):
        A = operators[i].get_dagger()
        for j in range(i+1,n):
            B = operators[j].get_dagger()
            for k in range(n):
                C = operators[k]
                for l in range(k+1,n):
                    D = operators[l]
                    Op = A*B*D*C
                    out = wf.dot(Op*wf) # overlap
                    cm[i,j,k,l] = out
#                    if np.abs(out)>1e-4:
#                       print(i,j,k,l,np.round(out,2))
#    print("Trace",np.sum([cm[i,i,i,i] for i in range(n)]))
    cm = four2two(cm) 
    if np.sum(np.abs(cm-np.conjugate(cm.T)))>1e-4: raise
    return cm # return matrix


def entropy_density(dm):
    """Compute the entropy density"""
    from scipy.linalg import eigh
    (es,ws) = eigh(dm) # compute
    ws = ws.T # transpose
    d = 0. # initialize
    for i in range(len(es)):
        e = es[i]
        w = ws[i] ; w2 = np.abs(w)**2
        if e>1e-8:
            d = d + e*np.log(e)*w2
    return -d




def four2two(m):
    n = m.shape[0] # dimension
    out = np.zeros((n*n,n*n),dtype=np.complex128) # output
    ii = 0
    for i in range(n-1):
      for j in range(i+1,n):
        jj = 0
        for k in range(n-1):
          for l in range(k+1,n):
              out[ii,jj] = m[i,j,k,l]
              jj += 1
        ii += 1
    return out[0:ii,0:ii]





def cpp_correlation_matrix(wf):
    """Compute the correlation matrix using the C++
    specialized function"""
    self = wf.MBO
    return self._session.correlation_matrix(wf.cpp_handle)


def get_four_correlation_tensor(wf,ctmode=None,**kwargs):
    """Return the correlation tensor as <Cdag_i C_j Cdag_k C_l>.

    ctmode=None (the default) auto-selects the fastest method actually
    available for this wavefunction's backend/chain type -- "sweep"
    whenever it applies (itensor_version in (3,"python"), non-native-
    spinful fermionic sites), else "full" whenever it applies, else the
    always-correct but slowest "explicit" fallback (see
    _four_correlation_tensor_default_ctmode()). Passing a ctmode
    explicitly is still a hard request: it raises rather than silently
    falling back if that method isn't available for this wavefunction."""
    if ctmode is None:
        ctmode = _four_correlation_tensor_default_ctmode(wf)
    if ctmode=="explicit":
        return get_four_correlation_tensor_explicit(wf,**kwargs)
    elif ctmode=="full":
        return get_four_correlation_tensor_cpp(wf,**kwargs)
    elif ctmode=="sweep":
        return get_four_correlation_tensor_sweep(wf,**kwargs)
    elif ctmode=="fold":
        return get_four_correlation_tensor_fold(wf,**kwargs)
    else: raise




def _single_factor_modes(ops):
    """[(operator_name, site_1based), ...] read off a list of one-factor
    MultiOperators (a chain's own `self.C`/`self.Cdag`), or None if any entry
    is not a single named operator on a single site.

    This is what lets `ctmode="fold"` be flavor-agnostic: a spinless chain's
    `C[m]` is `[[1.0, ["C", m]]]` and a native spinful chain's is
    `[[1.0, ["Cup", m//2]]]` / `[[1.0, ["Cdn", m//2]]]`, so the mode-to-
    (name, site) map falls straight out with no per-class special-casing."""
    modes = []
    for op in ops:
        terms = getattr(op, "op", None)
        if not terms or len(terms) != 1:
            return None
        term = terms[0]
        if len(term) != 2 or abs(complex(term[0]) - 1.0) > 1e-12:
            return None
        name, site = term[1]
        modes.append((name, int(site) + 1))
    return modes


def get_four_correlation_tensor_fold(wf, accelerate=True, **kwargs):
    """Four-point tensor by local operator folds
    (pyitensor.chain.Chain.four_correlation_tensor_fold).

    Works for any chain whose `C`/`Cdag` are single named operators on a
    single site -- spinless *and* native spinful -- and needs no MPO per
    tuple, unlike `ctmode="explicit"`, which is otherwise the only option
    for native spinful sites."""
    MBO = wf.MBO
    version = getattr(MBO, "itensor_version", None)
    if version not in ("python", 3):
        raise ValueError("ctmode='fold' needs itensor_version='python' or 3")
    cdag_ops = _single_factor_modes(MBO.Cdag)
    c_ops = _single_factor_modes(MBO.C)
    if cdag_ops is None or c_ops is None:
        raise ValueError("ctmode='fold': this chain's C/Cdag are not "
                          "single-site single-operator MultiOperators")
    session = MBO._session
    if version == 3:
        # the C++ binding takes parallel name/site vectors rather than a
        # list of pairs, which pybind11 converts without a custom caster
        return session.four_correlation_tensor_fold(
            wf.cpp_handle,
            [nm for nm, _s in cdag_ops], [s for _nm, s in cdag_ops],
            [nm for nm, _s in c_ops], [s for _nm, s in c_ops],
            accelerate)
    return session.four_correlation_tensor_fold(
        wf.cpp_handle, cdag_ops, c_ops, accelerate)


def get_four_correlation_tensor_explicit(wf,accelerate=True,**kwargs):
    """Compute the four field correlation tensor explicitly.

    <Cdag_i C_j Cdag_k C_l>^dagger = Cdag_l C_k Cdag_j C_i, so
    ct[i,j,k,l] and ct[l,k,j,i] are complex conjugates of each other --
    accelerate=True (the default) exploits this to only ever evaluate
    one representative of each (current,conjugate) pair, same as
    get_four_correlation_tensor_cpp()'s own accelerate flag and
    pyitensor.chain.Chain.four_correlation_tensor()'s (this was the one
    of the three implementations that didn't have it yet).
    """
#    C = [wf.MBO.toMPO(Ci) for Ci in wf.MBO.C] # operators
#    Cdag = [wf.MBO.toMPO(Cdagi) for Cdagi in wf.MBO.Cdag] # operators
    C = wf.MBO.C
    Cdag = wf.MBO.Cdag
    n = len(C) # dimension
    ct = np.zeros((n,n,n,n),dtype=np.complex128)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    current,conjugate = (i,j,k,l),(l,k,j,i)
                    if accelerate and current>conjugate: continue
                    Op = (Cdag[i]*C[j])*(Cdag[k]*C[l])
                    out = wf.aMb(Op,wf) # overlap
                    ct[i,j,k,l] = out
                    if accelerate and current!=conjugate:
                        ct[l,k,j,i] = np.conjugate(out)
    return ct


def get_four_correlation_tensor_cpp(wf,accelerate=True,**kwargs):
    """Compute the correlation tensor using the C++
    specialized function"""
    self = wf.MBO
    from .. import fermionchain
    if type(self)==fermionchain.Spinful_Fermionic_Chain_Native:
        # Native (Electron/Hubbard) sites: ITensor's ElectronSite has no
        # plain "Cdag"/"C" operator (only Cup/Cdn/Cdagup/Cdagdn), so the
        # generic four_correlation_tensor() C++ method can't be used
        # as-is -- four_correlation_tensor_spinful() is the flavor-
        # resolved counterpart (mpscpp3/chain_session.h), returning the
        # same tensor over the 2*n flat fermionic modes this chain's n
        # native sites represent.
        return self._session.four_correlation_tensor_spinful(wf.cpp_handle,accelerate)
    return self._session.four_correlation_tensor(wf.cpp_handle,accelerate)


def get_four_correlation_tensor_sweep(wf,accelerate=True,**kwargs):
    """Compute the four-point correlator tensor with the fast,
    single-sweep, environment-reuse method (mpscpp3/chain_session.h's
    Chain::four_correlation_tensor_sweep for itensor_version=3,
    pyitensor/chain.py's own port of the same algorithm for
    itensor_version="python" -- see either docstring for the algorithm,
    which follows the reuse idea of ITensorCorrelators.jl,
    https://github.com/ITensor/ITensorCorrelators.jl, rather than
    get_four_correlation_tensor_cpp()'s independent per-tuple AutoMPO
    build). Real, measured speedup over ctmode="full" that grows with
    system size (roughly 1.2x at n=6 up to >2.5x at n=16 for a Hubbard
    chain at maxm=60 under itensor_version=3 -- reproduced directly by
    examples/staticcorrelators/four_correlation_tensor_sweep_VS_full),
    at machine-precision agreement -- this is the default ctmode whenever
    it's available (see get_four_correlation_tensor()'s ctmode=None
    handling).

    accelerate only gates the repeated-index entries here,
    unlike get_four_correlation_tensor_cpp()'s own accelerate, which
    skips ~half of the *dominant* per-tuple AutoMPO builds via conjugate-
    pair symmetry -- there's no equivalent saving available in this
    method's dominant pairwise-distinct sweep (see either
    four_correlation_tensor_sweep()'s own docstring, C++ or pyitensor,
    for why), so don't expect the usual ~2x win from accelerate=True here.

    Only implemented for itensor_version in (3,"python") on the plain
    (non-native-spinful) fermionic sites: itensor_version=3 needs the
    compiled mpscpp3 extension, and neither has a counterpart yet for
    Spinful_Fermionic_Chain_Native (whose four_correlation_tensor_spinful()
    instead hands ITensor's AutoMPO the literal flavor-resolved
    Cdagup/Cup/Cdagdn/Cdn names and lets its own automatic fermionic
    threading handle two flavors sharing one physical site -- a different
    JW convention this sweep's per-flat-mode gap rule doesn't reproduce
    as-is)."""
    self = wf.MBO
    from .. import fermionchain
    if type(self)==fermionchain.Spinful_Fermionic_Chain_Native:
        raise NotImplementedError("ctmode=\"sweep\" has no native-spinful-site "
                "counterpart yet; use ctmode=\"full\" "
                "(four_correlation_tensor_spinful) instead")
    if self.itensor_version not in (3,"python") or \
            not hasattr(self._session,"four_correlation_tensor_sweep"):
        raise NotImplementedError("ctmode=\"sweep\" requires itensor_version=3 "
                "(with the compiled mpscpp3 extension) or \"python\", got "
                "itensor_version="+str(self.itensor_version))
    return self._session.four_correlation_tensor_sweep(wf.cpp_handle,accelerate)


def _four_correlation_tensor_default_ctmode(wf):
    """Pick the best available ctmode for get_four_correlation_tensor()
    when the caller didn't request one explicitly: "sweep" whenever it
    applies (itensor_version in (3,"python"), non-native-spinful); for
    native-spinful, "fold" under itensor_version="python" and "full" under
    3; else "full" whenever it applies (itensor_version in (2,3,"python"));
    else "explicit" (always correct, backend-agnostic, and by far the
    slowest -- it is a per-tuple MPO build and full-chain sweep). Only ever reached
    for DMRG-backed wavefunctions (mps.py/mpsjulialive's own
    get_four_correlation_tensor): ED-backed ones (edtk/edchain.py,
    pyfermion/mbfermion.py) always force ctmode="explicit" themselves and
    never call this -- so wf.MBO here is guaranteed a proper
    Many_Body_Chain with .itensor_version/._session, but getattr() is used
    anyway rather than assuming that, since a missing/unexpected attribute
    should just fall back to "explicit" (always correct) rather than
    crash."""
    self = wf.MBO
    from .. import fermionchain
    is_native_spinful = type(self)==fermionchain.Spinful_Fermionic_Chain_Native
    itensor_version = getattr(self,"itensor_version",None)
    session = getattr(self,"_session",None)
    if not is_native_spinful and itensor_version in (3,"python") \
            and hasattr(session,"four_correlation_tensor_sweep"):
        return "sweep"
    if is_native_spinful:
        # "fold" before "full": for native spinful sites the only previous
        # option under itensor_version="python" was "explicit", which builds
        # an MPO and sweeps the whole chain per tuple (profiled at 4 sites:
        # 55% to_mpo, 38% inner). Local folds are exact to machine precision
        # against it and measured 8-12x faster.
        if itensor_version in ("python", 3) \
                and hasattr(session,"four_correlation_tensor_fold"):
            return "fold"
        if itensor_version==3 and hasattr(session,"four_correlation_tensor_spinful"):
            return "full"
        return "explicit"
    if itensor_version in (2,3,"python") and hasattr(session,"four_correlation_tensor"):
        return "full"
    return "explicit"
