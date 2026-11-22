"""A single-particle basis in which a correlated state is weakly entangled,
from the Fishman-White circuit of a reference correlation matrix (Wu,
Kloss, Krinitsin, Fishman, Pixley and Stoudenmire, arXiv:2212.09798).

The circuit of circuit.correlation_matrix_to_circuit brings a Slater
determinant to a product of occupations. Its unitary V, the product of
the Givens rotations, defines new modes f_r^dag = sum_i V[r,i] c_i^dag,
so that c_i = sum_r V[r,i] f_r. Written in the f modes, the reference
state is a product state, and a correlated state close to it is weakly
entangled. For an impurity model the reference is the ground state
without the interaction, the gates touching the impurity are dropped
(the paper's choice), so that V is the identity on the impurity and the
interaction stays local, and the circuit runs from the far end of each
bath toward the impurity, so that the modes the circuit cannot decouple
end next to it.

V is unitary whatever the reference, so the rotated Hamiltonian has the
spectrum of the original one exactly: the reference only decides how
weakly entangled the ground state is in the new basis."""

import numpy as np

from . import circuit


def runs(n, fixed):
    """The maximal runs [a, b) of the modes 0..n-1 that are not fixed."""
    fixed = set(fixed)
    out = []
    a = None
    for i in range(n+1):
        free = i < n and i not in fixed
        if free and a is None:
            a = i
        elif not free and a is not None:
            out.append((a, i))
            a = None
    return out


def rotation(L, fixed=(), eigval_cutoff=1e-8, maxblocksize=None):
    """(V, info) for the correlation matrix L of one species of modes.

    V is unitary and the identity on the fixed modes. Each run of modes
    that are not fixed is rotated on its own, with the circuit of L on
    the run and the fixed mode next to it, started at the far end and
    with the gates on the fixed mode dropped. A run between two fixed
    modes has no far end, and raises. info is one dict per run: the
    modes, the number of gates kept and the largest distance to 0 or 1
    of the occupations the circuit decoupled."""
    L = np.array(L, dtype=complex)
    n = L.shape[0]
    fixed = sorted(set(int(f) for f in fixed))
    for f in fixed:
        if not 0 <= f < n:
            raise ValueError("fixed mode %d is out of range 0..%d" % (f, n-1))
    V = np.eye(n, dtype=complex)
    info = []
    for (a, b) in runs(n, fixed):
        left, right = a == 0, b == n
        if left and right:  # nothing fixed
            idx = list(range(n))
        elif left:  # the circuit runs from mode 0 toward the fixed mode b
            idx = list(range(a, b))+[b]
        elif right:  # from mode n-1 toward the fixed mode a-1
            idx = list(range(b-1, a-1, -1))+[a-1]
        else:
            raise NotImplementedError(
                "the modes %d to %d sit between two fixed modes, and a "
                "circuit that runs toward the impurity needs a far end; "
                "order the chain so that every bath is outside the fixed "
                "modes" % (a, b-1))
        stop = None if (left and right) else len(idx)-1  # the fixed mode
        ns, rots, err = circuit.correlation_matrix_to_circuit(
                L[np.ix_(idx, idx)], eigval_cutoff=eigval_cutoff,
                maxblocksize=maxblocksize)
        W = np.eye(len(idx), dtype=complex)
        kept = 0
        for (k, c, s) in rots:
            if stop is not None and stop in (k, k+1):
                continue  # keeps the interaction on the fixed mode local
            g = np.array([[c, s], [-np.conj(s), c]])
            W[k:k+2, :] = g @ W[k:k+2, :]
            kept += 1
        V[np.ix_(idx, idx)] = W
        info.append({"modes": (a, b), "gates": kept, "error": err})
    return V, info


# Operators allowed on a rotated mode, and the matrix that maps the
# coefficients of a factor of that kind onto the new modes: c_i =
# sum_r V[r,i] f_r and c_i^dag = sum_r conj(V[r,i]) f_r^dag.
_LADDER = {"C": False, "Cdag": True}


def rotate_terms(terms, V, rotated, scale_rule=1e-12):
    """The terms [c, [name, site], ...] of an operator written in the
    original modes, rewritten in the modes of V.

    A factor on a rotated mode must be C, Cdag or N (taken as Cdag C), at
    most two of them per term, since a quartic term on rotated modes
    becomes a sum over four indices. Terms that share every other factor
    and the slots of their rotated factors are collected into a vector
    or a matrix of coefficients, which is rotated in one product. A new
    coefficient below scale_rule times the largest one of its group is
    dropped, the threshold the MultiOperator filter applies anyway."""
    rotated = set(rotated)
    Vc = V.conj()
    n = V.shape[0]
    groups = {}
    out = []
    for term in terms:
        c = term[0]
        factors = []
        for (name, site) in term[1:]:
            if site in rotated and name == "N":
                factors += [["Cdag", site], ["C", site]]
            else:
                factors.append([name, site])
        slots = []
        for p, (name, site) in enumerate(factors):
            if site not in rotated or name == "Id":
                continue
            if name not in _LADDER:
                raise ValueError(
                    "the operator %s on site %d cannot be rotated: on a "
                    "rotated mode only C, Cdag and N can appear (make the "
                    "site fixed, or write the term with C and Cdag)"
                    % (name, site))
            slots.append(p)
        if not slots:
            out.append([c]+factors)
            continue
        if len(slots) > 2:
            raise ValueError(
                "a term with %d fermion operators on the rotated modes %s: "
                "an interaction must act on fixed modes only, since the "
                "rotation spreads it over the whole bath"
                % (len(slots), sorted({factors[p][1] for p in slots})))
        key = tuple((name, None) if p in slots else (name, site)
                    for p, (name, site) in enumerate(factors))
        sites = tuple(factors[p][1] for p in slots)
        if key not in groups:
            groups[key] = (slots, np.zeros((n,)*len(slots), dtype=complex))
        groups[key][1][sites] += c
    for key, (slots, M) in groups.items():
        maps = [Vc if _LADDER[key[p][0]] else V for p in slots]
        if len(slots) == 1:
            M = maps[0] @ M
        else:
            M = maps[0] @ M @ maps[1].T
        tol = scale_rule*np.max(np.abs(M))
        same = len(slots) == 2 and key[slots[0]][0] == key[slots[1]][0]
        for r in zip(*np.nonzero(np.abs(M) > tol)):
            if same and r[0] == r[1]:
                continue  # c_r c_r = 0
            factors = [list(f) for f in key]
            for p, site in zip(slots, r):
                factors[p][1] = int(site)
            coeff = M[r]
            if abs(coeff.imag) == 0.0:
                coeff = float(coeff.real)
            out.append([coeff]+factors)
    return out
