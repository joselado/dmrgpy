"""The Gaussian basis of an impurity model: bath orbitals rotated by the
Fishman-White circuit of a reference state, so that the correlated ground
state is weakly entangled (Wu, Kloss, Krinitsin, Fishman, Pixley and
Stoudenmire, arXiv:2212.09798).

The reference is a correlation matrix Lambda[i,j] = <c_i^dag c_j>, by
default the Slater determinant of a hopping matrix h with nf fermions,
for an impurity model the ground state without the interaction (the
resonant level model). Its circuit (gaussianmpstk/basis.py) gives a
unitary V and new modes f_r^dag = sum_i V[r,i] c_i^dag, written on the
same sites of the chain. The gates touching the fixed modes (the
impurity) are dropped, so V is the identity there and an interaction on
them stays as it is, and each bath is rotated from its far end toward
the impurity. The reference is then a product state except next to the
impurity, and the interacting ground state is close to it.

GaussianBasis.rotate rewrites an operator of the original modes in the
new ones, a Hamiltonian as well as an observable. A bilinear of the bath
becomes a sum over all pairs of new modes, so the rotated Hamiltonian
has longer-range hoppings and a larger MPO; the gain is a ground state
that needs a much smaller bond dimension.

GaussianBasis.natural_basis takes the reference from a solved state
instead: the correlation matrix of the ground state of the rotated chain,
taken back to the original modes, so that the circuit is that of the
natural orbitals of the interacting state. The paper leaves this open.
Measured on a 32-orbital Anderson impurity (v3, 12 sweeps, source state
at maxm=32, docs/gaussian_basis_comparison/README.md), against the
reference without U at the level ed+U/2: with eigval_cutoff=1e-5 the
error at maxm=8 is smaller in all four cases (uniform bath and Wilson
chain, at and away from particle-hole symmetry), by 1.4 to 6, and at
maxm=16 in three of the four; at maxm=32 it is about 2 times larger on
the uniform bath. The MPO is 26 to 46 against 26 to 34. With the 1e-8 of
a pure reference the uniform bath gains more (6 to 9 at maxm=8), but the
Wilson chain away from particle-hole symmetry loses 2.4 at maxm=8 (with
30 sweeps; 5 with 12, which do not converge that basis), with an MPO of
66. A second step gives nothing more. So it is a method to call,
not a default of any route.

On a Spinful_Fermionic_Chain, h, cm and fixed are given per orbital (an
n x n matrix is taken for both spins, as in get_gaussian_state, and
fixed lists orbitals, as fc.impurity does), and when the reference does
not mix the spins each spin is rotated on its own, so that the rotated
Hamiltonian keeps every conservation law of the original one."""

import numpy as np

from .gaussianmpstk import basis as gbasis

# Distance to [0, 1] above which an eigenvalue is not an occupation.
OCCUPATION_TOLERANCE = 1e-8

# eigval_cutoff of natural_basis: a measured state is not Gaussian, and
# with the 1e-8 of a pure reference the blocks of the circuit grow over
# the correlated part of the bath, which doubled the MPO of a Wilson chain
# and lost accuracy there (module docstring).
NATURAL_EIGVAL_CUTOFF = 1e-5


def _check_chain(chain):
    from .fermionchain import Fermionic_Chain, Majorana_Chain
    if not isinstance(chain, Fermionic_Chain) or isinstance(chain,
                                                            Majorana_Chain):
        raise TypeError("get_gaussian_basis needs a Fermionic_Chain or a "
                        "Spinful_Fermionic_Chain (spinless sites); the "
                        "native spinful chain is not supported")


def _fixed_modes(chain, fixed):
    from .gaussianmps import _is_spinful
    fixed = [int(f) for f in np.atleast_1d(fixed)] if fixed is not None \
        else []
    nmax = chain.ns//2 if _is_spinful(chain) else chain.ns
    for f in fixed:
        if not 0 <= f < nmax:
            raise ValueError("get_gaussian_basis: fixed %s %d is out of "
                             "range 0..%d" % ("orbital" if
                             _is_spinful(chain) else "mode", f, nmax-1))
    if _is_spinful(chain):
        return fixed, sorted([2*f for f in fixed]+[2*f+1 for f in fixed])
    return fixed, sorted(fixed)


def _reference(chain, h, nf, cm):
    """The reference correlation matrix on the chain's modes."""
    from . import gaussianmps
    if (h is None) == (cm is None):
        raise ValueError("get_gaussian_basis: give either h (with nf) or cm")
    if cm is None:
        return gaussianmps.correlation_matrix(chain, h, nf)
    if nf is not None:
        raise ValueError("get_gaussian_basis: nf goes with h, not cm")
    cm = gaussianmps._square(cm, "cm")
    if gaussianmps._is_spinful(chain) and 2*cm.shape[0] == chain.ns:
        cm = np.kron(cm, np.eye(2))
    if cm.shape[0] != chain.ns:
        raise ValueError("get_gaussian_basis: cm is %dx%d but the chain has "
                         "%d fermionic modes" % (cm.shape[0], cm.shape[1],
                                                 chain.ns))
    if np.max(np.abs(cm-cm.conj().T)) > 1e-12*gaussianmps._scale(cm):
        raise ValueError("get_gaussian_basis: cm must be Hermitian")
    ev = np.linalg.eigvalsh(cm)
    if ev[0] < -OCCUPATION_TOLERANCE or ev[-1] > 1.0+OCCUPATION_TOLERANCE:
        raise ValueError("get_gaussian_basis: cm has an eigenvalue %.3e "
                         "outside [0, 1], so it is not a correlation matrix "
                         "<Cdag_i C_j>" % (ev[0] if ev[0] < 0 else ev[-1]))
    return cm


def _rotation(chain, L, fixed_orbitals, fixed_modes, eigval_cutoff,
              maxblocksize):
    from .gaussianmps import _is_spinful
    kw = dict(eigval_cutoff=eigval_cutoff, maxblocksize=maxblocksize)
    if not _is_spinful(chain) or \
            np.max(np.abs(L[0::2, 1::2])) > 1e-12:  # the spins mix
        return gbasis.rotation(L, fixed_modes, **kw)
    Lup, Ldn = L[0::2, 0::2], L[1::2, 1::2]
    Vup, iup = gbasis.rotation(Lup, fixed_orbitals, **kw)
    if np.max(np.abs(Lup-Ldn)) <= 1e-12:
        Vdn, idn = Vup, iup
    else:
        Vdn, idn = gbasis.rotation(Ldn, fixed_orbitals, **kw)
    V = np.zeros_like(L, dtype=complex)
    V[0::2, 0::2], V[1::2, 1::2] = Vup, Vdn
    info = [dict(d, spin="up") for d in iup]+[dict(d, spin="dn") for d in idn]
    return V, info


class GaussianBasis:
    """The rotated modes f_r^dag = sum_i V[r,i] c_i^dag of a chain.

    V: the unitary, on the chain's modes.
    fixed: the fixed orbitals (spinful chain) or modes, as given.
    occupations: <f_r^dag f_r> of the reference, per mode; the modes far
        from 0 and 1 are the ones that carry the entanglement.
    info: one dict per rotated run of modes (and spin): the modes, the
        number of gates, and the largest distance to 0 or 1 of a
        decoupled occupation.
    spin_mixing: for a basis from natural_basis, the largest element of
        the spin-mixing block of the measured correlation matrix that was
        dropped to keep the spins apart (0 otherwise)."""

    def __init__(self, chain, V, fixed, fixed_modes, reference, info,
                 eigval_cutoff=1e-8, maxblocksize=None):
        self.chain = chain
        self.V = V
        self.fixed = fixed
        self.fixed_modes = fixed_modes
        self.info = info
        self.eigval_cutoff = eigval_cutoff
        self.maxblocksize = maxblocksize
        self.spin_mixing = 0.0
        self._reference = reference
        self._rotated = [i for i in range(chain.ns) if i not in
                         set(fixed_modes)]

    def correlation_matrix(self):
        """The reference correlation matrix <f_r^dag f_s> in the new
        modes, V Lambda V^dagger."""
        return self.V @ self._reference @ self.V.conj().T

    @property
    def occupations(self):
        return np.real(np.diag(self.correlation_matrix()))

    def mixes_spins(self):
        """Whether V mixes the two spins of a spinful chain."""
        from .gaussianmps import _is_spinful
        if not _is_spinful(self.chain):
            return False
        return bool(np.any(self.V[0::2, 1::2] != 0) or
                    np.any(self.V[1::2, 0::2] != 0))

    def original_correlation_matrix(self, wf=None, **kwargs):
        """<c_i^dag c_j> on the original modes of a state of the chain
        written in this basis (by default its ground state): with
        Lambda'[r,s] = <f_r^dag f_s> measured on the chain, it is
        V^dagger Lambda' V. Extra keywords go to get_correlation_matrix
        (mode="ED", for instance)."""
        cm = np.array(self.chain.get_correlation_matrix(wf=wf, **kwargs),
                      dtype=complex)
        return self.V.conj().T @ cm @ self.V

    def natural_basis(self, wf=None, eigval_cutoff=None, maxblocksize=None,
                      **kwargs):
        """The Gaussian basis of the correlation matrix of a state of the
        chain written in this basis (by default its ground state), with
        the same fixed orbitals: the circuit of the natural orbitals of
        the interacting state instead of those of the reference. The new
        basis is relative to the original modes, so it rotates the
        original Hamiltonian, as this one does.

        When this basis keeps the spins apart, the spin-mixing block of
        the measured matrix is dropped, so that the new one keeps them
        apart too (its largest element is kept as spin_mixing).
        eigval_cutoff defaults to NATURAL_EIGVAL_CUTOFF and maxblocksize
        to that of this basis; extra keywords go to
        get_correlation_matrix."""
        cm = self.original_correlation_matrix(wf=wf, **kwargs)
        cm = 0.5*(cm+cm.conj().T)
        dropped = 0.0
        if _is_spinful_chain(self.chain) and not self.mixes_spins():
            dropped = float(max(np.max(np.abs(cm[0::2, 1::2])),
                                np.max(np.abs(cm[1::2, 0::2]))))
            cm[0::2, 1::2] = 0.0
            cm[1::2, 0::2] = 0.0
        if eigval_cutoff is None:
            eigval_cutoff = NATURAL_EIGVAL_CUTOFF
        if maxblocksize is None:
            maxblocksize = self.maxblocksize
        b = gaussian_basis(self.chain, cm=cm, fixed=self.fixed,
                           eigval_cutoff=eigval_cutoff,
                           maxblocksize=maxblocksize)
        b.spin_mixing = dropped
        return b

    def rotate(self, op):
        """The operator op of the original modes, written in the new ones:
        a Hamiltonian to set with set_hamiltonian, or an observable of the
        original bath to measure in the new basis. On a rotated mode only
        C, Cdag and N may appear, at most two per term; an interaction
        must act on the fixed modes only."""
        from . import multioperator
        op = multioperator.obj2MO(op)
        terms = gbasis.rotate_terms(multioperator._filter_small(op.op),
                                    self.V, self._rotated,
                                    scale_rule=multioperator.clean_threshold)
        out = multioperator.MultiOperator(term=False)
        out.op = terms
        out.i = len(terms)-1
        return out


def gaussian_basis(chain, h=None, nf=None, cm=None, fixed=None,
                   eigval_cutoff=1e-8, maxblocksize=None):
    """The Gaussian basis of the reference state given by the hopping
    matrix h and the number of fermions nf (its Slater determinant), or
    by a correlation matrix cm (<Cdag_i C_j>, eigenvalues in [0, 1], not
    necessarily a pure state), with the fixed orbitals (spinful chain) or
    modes left out of the rotation.

    eigval_cutoff: how close to 0 or 1 an occupation must be for its mode
        to be decoupled; it sets how close to a product state the
        reference is in the new basis.
    maxblocksize: the largest block of the circuit (default: no limit).
        A smaller one gives fewer long-range hoppings in the rotated
        Hamiltonian and a reference that is less disentangled."""
    _check_chain(chain)
    fixed, fixed_modes = _fixed_modes(chain, fixed)
    L = _reference(chain, h, nf, cm)
    V, info = _rotation(chain, L, fixed, fixed_modes, eigval_cutoff,
                        maxblocksize)
    return GaussianBasis(chain, V, fixed, fixed_modes, L, info,
                         eigval_cutoff=eigval_cutoff,
                         maxblocksize=maxblocksize)


def _is_spinful_chain(chain):
    from .gaussianmps import _is_spinful
    return _is_spinful(chain)
