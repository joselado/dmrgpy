"""Exact diagonalization of a chain that mixes native spin sites and
spinless fermion modes: the ED backend of
mixedchain.Mixed_Spin_Fermion_Chain.

The Hilbert space is the tensor product of the physical sites in chain
order. A fermion mode (type code 0) is a two-state site, |0> and |1>; a
spin S site (type code 2S+1, as in spinchain.label2site) has the 2S+1
states m = S, S-1, ..., -S. A fermion operator carries a Jordan-Wigner
string over the fermion modes to its left, and only over those:

    C_p = prod_{k<p, k a fermion mode} (1 - 2 n_k) a_p ,

so a spin site between two fermion modes gives no sign. This is the
convention of the DMRG backends: multioperatortk/jordanwigner.py puts the
string F on every site to the left of p, and F is the identity on a spin
site (ITensor v3's SiteSet::op() fallback, and pyitensor). With it the C_p
obey the canonical anticommutation relations among themselves and commute
with every spin operator, which is what a spin next to a conduction
electron means physically.

Operators are named per physical site, with the names the DMRG site types
use: C, Cdag, N, A, Adag (no string), F and Id on a fermion mode; Sx, Sy,
Sz, Sp/S+, Sm/S-, S2, Sx2, Sy2, Sz2, F (the identity) and Id on a spin
site. Everything else (ground state, vev, excited states, dynamical
correlators, time evolution, conserved sectors) is the generic EDchain
machinery acting on these matrices.

Checked in tests/test_mixed_chain_ed.py: an all-spin chain gives exactly
the spectrum of spinchain.Spin_Chain's ED, an all-fermion chain exactly
the spectrum of fermionchain.Spinful_Fermionic_Chain's ED, and genuinely
mixed chains (a spin 1/2 or a spin 1 between fermion locations, with a
hopping across the spin) agree with the independent dense oracle of
tests/_dense_jw_reference.py and with DMRG on itensor_version=3 and
"python".
"""
import numpy as np
import scipy.sparse as sparse
from .edchain import EDchain
from .. import multioperator

maxdim = 2**20 # largest Hilbert space this backend builds


def spin_matrices(dim):
    """Sx, Sy, Sz, S+ and S- of a spin S=(dim-1)/2, in the basis
    m = S, S-1, ..., -S"""
    S = (dim-1)/2.
    ms = S - np.arange(dim) # m on each basis state
    sp = np.zeros((dim,dim),dtype=np.complex128)
    for k in range(1,dim): # <m+1|S+|m>, m = ms[k], m+1 = ms[k-1]
        sp[k-1,k] = np.sqrt(S*(S+1) - ms[k]*(ms[k]+1))
    sm = sp.T.copy()
    sx = 0.5*(sp + sm)
    sy = -0.5j*(sp - sm)
    sz = np.diag(ms).astype(np.complex128)
    return sx,sy,sz,sp,sm


def _local_dimension(code):
    """Local dimension of a physical site from its type code"""
    if code==0: return 2 # spinless fermion mode
    if isinstance(code,(int,np.integer)) and code>=2: return int(code) # 2S+1
    raise ValueError("the mixed-chain ED backend builds fermion modes (type "
            "code 0) and spins (type code 2S+1 >= 2), got type code %r"
            % (code,))


def _kron(mats):
    """Sparse Kronecker product of a list of local matrices"""
    out = sparse.csc_matrix(np.array([[1.0]],dtype=np.complex128))
    for m in mats:
        out = sparse.kron(out,sparse.csc_matrix(m,dtype=np.complex128),
                format="csc")
    out.eliminate_zeros()
    return out


class MixedSpinFermionED(EDchain):
    """ED object for a list of physical site type codes, 0 for a spinless
    fermion mode and 2S+1 for a spin S (the same list the DMRG session of
    a Mixed_Spin_Fermion_Chain is built from)"""
    def __init__(self,site_types):
        super().__init__()
        self.site_types = list(site_types)
        self.ns = len(self.site_types)
        self.localdim = [_local_dimension(c) for c in self.site_types]
        self.fermionic = [c==0 for c in self.site_types]
        dim = int(np.prod(self.localdim))
        if dim>maxdim:
            raise ValueError("Hilbert space dimension %d exceeds the maximum "
                    "allowed for ED (%d); use DMRG instead" % (dim,maxdim))
        self.size = dim
        self.hamiltonian = None # MultiOperator, set by the chain
        self.wf0 = None # ground state
        self.e0 = None # ground state energy
        self._ids = [np.identity(d,dtype=np.complex128) for d in self.localdim]
        # local matrices of each site type
        a = np.array([[0.,1.],[0.,0.]],dtype=np.complex128) # a|1> = |0>
        self._fermion_local = {"A":a, "Adag":a.T.copy(),
                "N":np.diag([0.,1.]).astype(np.complex128),
                "F":np.diag([1.,-1.]).astype(np.complex128)}
        self._spin_local = dict()
    def _spin_ops(self,dim):
        """Local spin matrices of a spin site of dimension dim"""
        if dim not in self._spin_local:
            sx,sy,sz,sp,sm = spin_matrices(dim)
            ops = {"Sx":sx,"Sy":sy,"Sz":sz,"Sp":sp,"Sm":sm,
                    "Sx2":sx@sx,"Sy2":sy@sy,"Sz2":sz@sz,
                    "S2":sx@sx+sy@sy+sz@sz,
                    "F":np.identity(dim,dtype=np.complex128)}
            ops["X"],ops["Y"],ops["Z"] = sx,sy,sz
            ops["S+"],ops["S-"] = sp,sm
            self._spin_local[dim] = ops
        return self._spin_local[dim]
    def _build(self,name,i):
        """Many-body matrix of the single-site name at physical site i"""
        if name=="Id": return self.get_identity()
        if not isinstance(i,(int,np.integer)) or not 0<=i<self.ns:
            raise ValueError("the mixed-chain ED backend has no physical "
                    "site %r (it has %d)" % (i,self.ns))
        mats = list(self._ids)
        if self.fermionic[i]:
            loc = self._fermion_local
            if name in ("N","density","n"): mats[i] = loc["N"]
            elif name in ("A","Adag","F"): mats[i] = loc[name]
            elif name in ("C","Cdag"):
                mats[i] = loc["A"] if name=="C" else loc["Adag"]
                for k in range(i): # Jordan-Wigner string, fermion modes only
                    if self.fermionic[k]: mats[k] = loc["F"]
            else:
                raise ValueError("the mixed-chain ED backend has no operator "
                        "%r on the fermion mode %d (it builds C, Cdag, N, A, "
                        "Adag, F and Id there)" % (name,i))
        else:
            ops = self._spin_ops(self.localdim[i])
            if name in ops: mats[i] = ops[name]
            else:
                raise ValueError("the mixed-chain ED backend has no operator "
                        "%r on the spin site %d (it builds Sx, Sy, Sz, Sp, "
                        "Sm, S2, Sx2, Sy2, Sz2, F and Id there)" % (name,i))
        return _kron(mats)
    def get_operator(self,name,i=0):
        """Return an operator: a MultiOperator or an EDOperator as a
        (sector-restricted) many-body matrix, or a single-site name at
        physical site i as a full-space matrix"""
        from .edchain import EDOperator
        if type(name)==multioperator.MultiOperator:
            return self.sector_restrict(multioperator.MO2matrix(name,self),name)
        elif type(name)==EDOperator: return name.SO
        elif isinstance(name,str):
            key = (name,i)
            if key not in self.operators: self.operators[key] = self._build(name,i)
            return self.operators[key]
        else:
            raise TypeError("the mixed-chain ED backend cannot build an "
                    "operator from %r (%s): expected a MultiOperator, an "
                    "EDOperator or a single-site name"
                    % (name,type(name).__name__))
