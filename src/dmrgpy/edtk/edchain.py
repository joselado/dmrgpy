from ..algebra import algebra
from .. import multioperator
import scipy.sparse as sparse
import scipy.sparse.linalg as slg
from .one2many import one2many
import numpy as np


def _operator_label(name):
    """A short, readable description of the operator a sector rejected.
    MultiOperator has no __repr__ of its own, and the default
    "<...object at 0x...>" says nothing, so name the terms instead."""
    if name is None: return "unnamed"
    if isinstance(name,str): return name
    op = getattr(name,"op",None)
    if op is None: return type(name).__name__
    terms = []
    for iop in op[0:4]: # a few terms is enough to recognize it by
        terms.append("*".join(str(t[0])+"["+str(t[1])+"]" for t in iop[1:]))
    out = " + ".join(t for t in terms if len(t)>0)
    if len(op)>4: out = out+" + ..."
    return out if len(out)>0 else type(name).__name__


def _diagonal_of(Q,name):
    """Diagonal of a conserved-charge operator, as a real array over the
    ED basis. A conserved quantity of a chain is diagonal in the
    occupation/Sz product basis by construction; anything else cannot
    define a sector as a set of basis states, so say so rather than
    quietly using its diagonal."""
    # never densified: the charge operator lives on the *full* Hilbert
    # space, so Q.todense() here would be a dim x dim complex array (2 GB
    # at the 14-site chains ED already handles comfortably) for something
    # that only ever needs its diagonal
    if sparse.issparse(Q):
        d = np.asarray(Q.diagonal()).ravel().copy()
        R = (Q - sparse.diags(d,format="csr")).tocoo()
        offdiag = np.max(np.abs(R.data)) if R.nnz>0 else 0.0
    else:
        Q = np.asarray(Q)
        d = np.diag(Q).copy()
        offdiag = np.max(np.abs(Q-np.diag(d)))
    if offdiag>1e-8:
        raise ValueError("the conserved quantity "+str(name)+" is not "
                "diagonal in the ED basis, so it does not define a sector")
    if np.max(np.abs(d.imag))>1e-8:
        raise ValueError("the conserved quantity "+str(name)+" is not real")
    return d.real


class EDchain():
    """Generic class for an ED chain"""
    def __init__(self):
        self.operators = dict() # empty dictionary
        self.localdim = [] # empty list
        self.computed_gs = False
        self.Identity = None # placeholder for identity
        self.EDHamiltonian = None # initialize as None
        self.Diagonalized_Hamiltonian = None # initialize as None
        self._kpm_emax_cache = None # (e0,emax) cache for the KPM dynamical correlator
        # Conserved-sector (quantum-number) targeting, off by default --
        # see set_conserved_sector() below and Many_Body_Chain.
        # set_conserved_sector(), which is the public entry point.
        self.conserved_sector = None # dict {name: target}, or None
        self._sector_mask = None # boolean array over the full ED basis
        self._sector_charge = None # dict {name: charge per basis state}
    ##########################################
    # conserved-sector (quantum number) mode #
    ##########################################
    def set_conserved_sector(self,qns,charges=None):
        """Confine this ED chain to one quantum-number sector.

        `qns` is a dict {name: target} (or None/{} to switch the sector
        off) and `charges` a dict {name: MultiOperator} giving the
        many-body operator each conserved quantity is measured by -- both
        supplied by the Many_Body_Chain this object was built for (see
        Many_Body_Chain.get_sector_charge_operators).

        The mechanism is the one the sector tests already use as their
        reference: a conserved charge is diagonal in the ED product
        basis, so a sector is a *set of basis states* and confining the
        calculation to it is taking a submatrix. Every many-body operator
        this object assembles from a MultiOperator is restricted to those
        basis states before it is handed out, so the eigensolver, the
        expectation values, the correlators and the dynamics all run
        inside the sector and nowhere else.

        Unlike the DMRG backends, the *full* Hilbert space is still built
        on the way there (the restriction happens after assembly), so this
        buys a smaller eigenproblem rather than a smaller construction.
        """
        # any cached matrix/spectrum was built for the previous sector
        self.EDHamiltonian = None
        self.Diagonalized_Hamiltonian = None
        self._kpm_emax_cache = None
        self.computed_gs = False
        self.conserved_sector = None # cleared first: the charge operators
        self._sector_mask = None     # below must be assembled unrestricted
        self._sector_charge = None
        if not qns: return # sector targeting off, nothing else to do
        charges = charges or dict()
        if not charges:
            raise NotImplementedError(
                "this ED backend defines no conserved quantities, so it "
                "cannot target the sector "+str(dict(qns)))
        unknown = [k for k in qns if k not in charges]
        if len(unknown)>0:
            raise ValueError(
                "unknown conserved quantity "+str(unknown)+" for this chain; "
                "it offers "+str(sorted(charges)))
        sel = None # boolean mask over the full basis
        cs = dict() # charge of every basis state, per quantity
        for k in sorted(qns):
            Q = multioperator.MO2matrix(charges[k],self) # full space
            q = _diagonal_of(Q,k) # values on the basis states
            cs[k] = q
            si = np.abs(q-float(qns[k]))<1e-6
            sel = si if sel is None else (sel & si)
        if not sel.any():
            raise ValueError(
                "empty sector "+str(dict(qns))+": no basis state of this "
                "chain carries those quantum numbers")
        self._sector_charge = cs
        self._sector_mask = sel
        self.conserved_sector = {k:qns[k] for k in qns}
    def sector_dimension(self):
        """Dimension of the space this chain actually works in -- the size
        of the conserved sector if one is set, the full Hilbert space
        otherwise"""
        if self._sector_mask is None: return self.get_identity().shape[0]
        return int(np.sum(self._sector_mask))
    def sector_nonconserving(self,A,tol=1e-10):
        """Names of the conserved quantities the many-body matrix A does
        not commute with (empty when no sector is set). The charges are
        diagonal, so this is just "does A connect two basis states of
        different charge", checked on A's own nonzero entries."""
        if self._sector_charge is None: return []
        M = A.tocoo() if sparse.issparse(A) else sparse.coo_matrix(np.asarray(A))
        if M.nnz==0: return []
        data = np.asarray(M.data)
        # relative threshold: a term that cancelled exactly against
        # another one (the S+S+/S-S- pair of Sx.Sx+Sy.Sy) leaves rounding
        # dust behind, and that dust is not a symmetry violation
        big = np.abs(data)>tol*max(np.max(np.abs(data)),1.0)
        rows,cols = np.asarray(M.row)[big],np.asarray(M.col)[big]
        out = []
        for k in sorted(self._sector_charge):
            q = self._sector_charge[k]
            if np.any(np.abs(q[rows]-q[cols])>1e-9): out.append(k)
        return out
    def sector_restrict(self,A,name=None):
        """Restrict an assembled many-body operator to the conserved
        sector, refusing the ones that do not conserve it.

        Refusing rather than silently projecting matches what the DMRG
        backends do with a flux-violating operator, and it is not merely
        cosmetic: P.A.P is exact for a static expectation value but
        identically *zero* for a charge-changing operator, so a dynamical
        correlator of C/S+ inside a fixed-charge sector would come back as
        a clean, wrong zero instead of an error."""
        if self._sector_mask is None: return A # no sector, nothing to do
        bad = self.sector_nonconserving(A)
        if len(bad)>0:
            raise ValueError(
                "this operator ("+_operator_label(name)+") does not conserve "
                +str(bad)+", which the conserved sector "
                +str(self.conserved_sector)+" fixes. Measure it after "
                "clearing the sector with set_conserved_sector().")
        sel = np.where(self._sector_mask)[0]
        return A[sel,:][:,sel] # the sector's submatrix
    def sector_embed(self,v):
        """Scatter a vector living in the conserved sector back into the
        full Hilbert space, padding with the zeros the sector forbids"""
        if self._sector_mask is None: return v
        out = np.zeros(self._sector_mask.shape[0],dtype=np.complex128)
        out[self._sector_mask] = v
        return out
    def get_operator(self,name,i=0):
        """Return an operator"""
        if type(name)==multioperator.MultiOperator: # input is a MO
            return self.sector_restrict(multioperator.MO2matrix(name,self),name)
        elif type(name)==str: # string
            return self.operators[(name,i)] # return the operator
        elif type(name)==EDOperator:
            return name.SO
        else: # unrecognized type
            print("Unrecognized operator in EDchain",type(name))
    def get_hamiltonian(self):
        """Return the Hamiltonian"""
        if self.EDHamiltonian is None:
            out = self.get_operator(self.hamiltonian) # return operator
            self.EDHamiltonian = out
            return out
        else: return self.EDHamiltonian
    def get_diagonalized_hamiltonian(self):
        """Return Hamiltonian in diagonal form, eigenvalues and eigenvectors"""
        if self.Diagonalized_Hamiltonian is None:
            h = self.get_hamiltonian() # get the Hamiltonian
            emu,vs = algebra.eigh(h) # eigenvectors and eigenvalues
            self.Diagonalized_Hamiltonian = (emu,vs) # store
            return (emu,vs)
        else: return self.Diagonalized_Hamiltonian
    def gs_energy(self):
        """Return ground state energy"""
        return self.get_excited(n=1)[0]
    def get_gs(self,array_mode=True):
        """Get the ground state"""
#        if array_mode: 
#            print("This will be deprecated")
#            return self.get_gs_array()
#        else: 
        return State(self.get_gs_array(),self)
    def get_gs_array(self):
        """Get ground state wavefunction"""
        if self.computed_gs: return self.wf0
        else: 
#          print("Computing GS")
#          e0,wf0 = algebra.ground_state(self.get_hamiltonian())
          (es,ws) = algebra.lowest_states(self.get_hamiltonian(),n=3)
          e0 = es[0]
          wf0 = ws[0]
          self.wf0 = wf0
          self.e0 = e0
          self.computed_gs = True
          return self.wf0
    def vev(self,op,T=0.,**kwargs):
        """Return a vacuum expectation value"""
        if T==0.: # zero temperature
            wf0 = self.get_gs_array()
            op = self.MO2matrix(op) # return operator
            return algebra.braket_wAw(wf0,op)
        else: # finite temperature
            from ..vevtk.thermalvev import thermal_vev_ex
            return thermal_vev_ex(self,op,T=T,**kwargs) # return thermal VEV
    def get_excited(self,**kwargs):
        """Excited states"""
        h = self.get_hamiltonian()
        return algebra.lowest_eigenvalues(h,**kwargs)
    def get_excited_states(self,**kwargs):
        """Excited states"""
        h = self.get_hamiltonian()
        (es,ws) = algebra.lowest_states(h,**kwargs)
        ws = [State(w,self) for w in ws] # transform to states
        return (es,ws)
    def create_operator(self,a,i=None,name=None):
        """Create the different operators"""
        ns = self.localdim # local dimensions
        if name is None: raise
        if i is None: raise
        if a.shape[0]!=ns[i]: raise
        # create operators in each site
        ids = [np.identity(n,dtype=np.complex128) for n in ns] # identities
        op = one2many(ids,a,i) # one to many body
        self.operators[(name,i)] = op # store in the dictionary
    def get_identity(self):
        if self.Identity is None:
            ids = [np.identity(n,dtype=np.complex128) for n in self.localdim] 
            op = one2many(ids,ids[0],0) # one to many body
            self.Identity = op
            return op
        else: return self.Identity
    def get_dynamical_correlator(self,**kwargs):
        from . import dynamics
        return dynamics.get_dynamical_correlator(self,**kwargs)
    def get_distribution(self,**kwargs):
        """Return a certain distribution"""
        from . import distribution
        return distribution.get_distribution(self,**kwargs)
    def exponential(self,h,wf):
        """Exponential of a wavefunction"""
#        wf = State(wf,self) # convert to State
        h = self.MO2matrix(h) # convert to matrix 
#        print(wf) ; exit()
        return State(algebra.expm(h)@wf.v,self) # return
    def get_ED_obj(self): return self
    def MO2matrix(self,m): 
        if type(m)==EDOperator: return m.SO # static operator
        else: return self.sector_restrict(multioperator.MO2matrix(m,self),m)
    def overlap(self,wf1,wf2):
        return wf1.dot(wf2) 
    def applyoperator(self,A,wf): 
        wf = State(wf,self)
        return A*wf #return self.MO2matrix(A)@wf
    def random_state(self):
        """Return a random state"""
        n = self.sector_dimension() # dimension (of the sector, if one is set)
        v = np.random.random(n)-.5 + 1j*(np.random.random(n)-.5)
        return State(v,self).normalize() # return the state
    def is_zero_operator(self,A):
        """Check if this is a zero operator"""
        return algebra.is_zero_matrix(self.obj2matrix(A))
    def obj2matrix(self,a):
        return obj2matrix(self,a)
    def applyinverse(self,A,wf):
        """Apply inverse, A is the matrix, wf the wavefunction"""
        return wf.applyinverse(A)




class State():
    """This is a dummy class to contain states"""
    def __init__(self,v,MBO):
        if type(v)==State:
            self.v = v.v.copy()
            self.MBO = v.MBO
        else:
            self.v = v # store the vector
            self.MBO = MBO # store the many-body object
        self.mode = "ED" # exact diagonalization mode
    def __rmul__(self,a):
        """Multiply by something"""
        from ..algebra.algebra import ismatrix
        if type(a)==multioperator.MultiOperator: # multioperator
            A = self.MBO.MO2matrix(a)  # get the matrix
            w = A@self.v # multiply
            return type(self)(w,self.MBO) # create a new object
        elif multioperator.isnumber(a):
            w = a*self.v # multiply
            return type(self)(w,self.MBO) # create a new object
        elif ismatrix(a):
            if a.shape[1]!=self.v.shape[0] \
                    and getattr(self.MBO,"conserved_sector",None):
                # a many-body matrix built on the full Hilbert space met a
                # wavefunction confined to a sector: the usual cause is an
                # operator that changes the conserved charge (C, S+), which
                # has no representation inside the sector at all
                raise ValueError(
                    "this matrix acts on the full Hilbert space (dimension "
                    +str(a.shape[1])+") but the wavefunction lives in the "
                    "conserved sector "+str(self.MBO.conserved_sector)
                    +" (dimension "+str(self.v.shape[0])+"). Clear the sector "
                    "with set_conserved_sector() to apply it.")
            w = a@self.v
            return type(self)(w,self.MBO) # create a new object
        else: raise # not implemented
    def __add__(self,a):
        if isinstance(a,State): return type(self)(self.v + a.v,self.MBO)
        else: raise
    def get_fermionic_parity(self,**kwargs):
        from ..fermionicparity import get_fermionic_parity
        return get_fermionic_parity(self,**kwargs) # parity of the state
    def get_conjugate(self):
        out = self.copy()
        out.v = np.conjugate(self.v) # conjugate wavefunction
        return out
    def __mul__(self,x):
        if multioperator.isnumber(x): 
            return type(self)(x*self.v,self.MBO)
        else: raise
    def __truediv__(self,a):
        if multioperator.isnumber(a): # number
          return (1./a)*self
        else: raise
    def __sub__(self,a):
        return self + (-1)*a
    def __neg__(self):
        return (-1)*self
    def overlap(self,a):
        if isinstance(a,State): # state object
            return np.dot(np.conjugate(self.v),a.v)
        else: raise # not implemented
    def copy(self):
        from copy import deepcopy
        return deepcopy(self)
    def dot(self,a):
        return np.sum(np.conjugate(self.v)*a.v)
    def aMb(self,M,b): return self.dot(M*b)
    def normalize(self):
        norm = np.sqrt(self.dot(self).real) # norm
        if norm>1e-8: return self/norm
        else: return None
    def get_correlation_entropy(self,**kwargs):
        from .. import entanglement
        return entanglement.get_correlation_entropy_from_wf(self,**kwargs)
    def get_four_correlation_tensor(self,**kwargs):
        """ED has no MPS to sweep or fold, so only ctmode="explicit"
        applies. It is forced rather than defaulted -- but via kwargs, so a
        caller who passes ctmode explicitly gets a clear error (or, for
        "explicit", just works) instead of TypeError: got multiple values
        for keyword argument 'ctmode', which is what a hardcoded keyword
        used to raise for any caller that named the mode at all."""
        from .. import entanglement
        ctmode = kwargs.pop("ctmode", "explicit")
        if ctmode not in (None, "explicit"):
            raise ValueError(
                "get_four_correlation_tensor: ctmode={!r} is not available "
                "for an ED-backed wavefunction; only \"explicit\" is".format(
                    ctmode))
        return entanglement.get_four_correlation_tensor(self,
                                    ctmode="explicit",**kwargs)
    def applyinverse(self,a,**kwargs):
        if type(a)==multioperator.MultiOperator: # multioperator
            A = self.MBO.MO2matrix(a)  # get the matrix
        elif type(a)==EDOperator: # multioperator
            A = a.SO  # get the matrix
        else: raise
        w = algebra.applyinverse(A,self.v)
        return type(self)(w,self.MBO) # create a new object



def obj2matrix(self,a):
    """ Transform to a matrix, self is an MBO and a the object to transform"""
    if type(a)==multioperator.MultiOperator: # multioperator
        return self.MO2matrix(a)  # get the matrix
    elif type(a)==EDOperator: # EDoperator
        return a.SO  # get the matrix
    else: raise




class EDOperator():
    """This is a dummy class for operators, so that it resembles the
    tensor network Static Operator. Useful for testing purposes"""
    def __init__(self,MO,MBO):
        """Init, takes as input a multioperator and the MBO"""
        self.MBO = MBO # store the many-body object
        if type(MO)==multioperator.MultiOperator: # multioperator
            self.SO = MBO.get_operator(MO) # generate the static operator
        elif type(MO)==EDOperator:
            self.SO = MO.SO.copy() # dummy copy
        else:
            # a bare `raise` here meant "RuntimeError: No active exception
            # to reraise" for the one thing that plausibly lands here: an
            # operator built for a different backend (a StaticOperator from
            # toMPO(mode="DMRG")) passed to a mode="ED" call
            raise TypeError(
                "EDOperator takes a MultiOperator or another EDOperator, "
                "got "+type(MO).__name__+". Operators built for the DMRG "
                "backend (toMPO(mode=\"DMRG\") -> StaticOperator) cannot "
                "be reused under mode=\"ED\"; build them with "
                "toMPO(mode=\"ED\"), or pass the MultiOperator itself.")
    def __add__(self,a):
        if multioperator.isnumber(a): # adding a number
            out = self.copy() # make a copy
            # identity of the operator's own dimension, which is the
            # conserved sector's rather than the full Hilbert space's
            # whenever this operator was built inside one
            o = self.SO + a*sparse.identity(self.SO.shape[0],
                                            dtype=np.complex128)
            out.SO = o.copy()
            return out # return 
        elif type(a)==EDOperator:
            out = self.copy() # make a copy
            out.SO = self.SO + a.SO
            return out
        else:
            print("Not recognized",type(a))
            raise
    def __radd__(self,a): return self + a # commutative
    def __sub__(self,a): return self + (-1)*a # commutative
    def __rsub__(self,a): return -self + a # commutative
    def __neg__(self): return (-1)*self  # use product
    def __mul__(self,v):
        from ..multioperator import MultiOperator
        if isinstance(v,State): # input is a state
            return State(self.SO@v.v,self.MBO)
        elif type(v)==EDOperator: # input is an MPO
            out = self.copy() # copy
            out.SO = self.SO@v.SO # matrix multiplication
            return out
        elif multioperator.isnumber(v): # adding a number
            out = self.copy() # copy
            out.SO = v*self.SO # scalar multiplication
            return out
        elif type(v)==MultiOperator: # input is a multioperator
            return self*EDOperator(v,self.MBO)
        else: 
            print("Incompatible object",type(v))
            raise
    def __rmul__(self,v):
        from ..multioperator import MultiOperator
        if type(v)==MultiOperator: # input is a multioperator
            return EDOperator(v,self.MBO)*self
        elif multioperator.isnumber(v): # adding a number
            return self*v
        else: raise
    def copy(self):
        from copy import deepcopy
        return deepcopy(self)
    def trace(self):
        from ..algebra.algebra import trace
        return trace(self.SO)
    def get_dagger(self):
        out = self.copy()
        from ..algebra.algebra import dagger
        out.SO = dagger(self.SO)
        return out









