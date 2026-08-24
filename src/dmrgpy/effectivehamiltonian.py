"""Effective low-energy Hamiltonians by projection onto a few eigenstates.

The construction is: take the `n` lowest eigenstates of `self.hamiltonian`,
project the Hamiltonian onto the manifold they span, and fit that
projected matrix as a real linear combination of the same projections of a
chosen set of operators. The fitted coefficients are the effective
couplings.

Three things about this are easy to get wrong, and are handled explicitly
here rather than left to the caller:

* **The manifold must be closed.** Cutting `n` inside a degenerate
  multiplet leaves a manifold that is not invariant under the model's own
  symmetries, and the fit then returns plausible-looking couplings that
  mean nothing (on a Hubbard dimer, whose low manifold is a singlet plus a
  3-fold triplet, `n=3` returned couplings ~50x the correct ones). Every
  entry point here checks the cut and raises rather than returning that.
* **The fit is linear**, so it is solved as a linear least-squares problem
  with a rank check, not with a nonlinear optimizer from a random start.
  See `fit_matrix`.
* **The eigenstates must be orthonormal.** ED's are to machine precision;
  DMRG's overlap-penalty excited states are not, so the basis is
  Lowdin-orthonormalized when it needs to be. See `orthonormalize`.
"""
import numpy as np


# Relative tolerance for calling two consecutive eigenvalues degenerate,
# i.e. for deciding that a manifold cut falls inside a multiplet. Scaled
# by the manifold's own energy spread, so it is dimensionless in practice.
DEGENERACY_TOL = 1e-6

# Above this deviation from the identity, the eigenstate Gram matrix is
# considered non-orthonormal and the basis is Lowdin-orthonormalized.
ORTHONORMALITY_TOL = 1e-10


def get_manifold(self, n, mode="DMRG", degeneracy_tol=DEGENERACY_TOL):
    """Return the `n` lowest eigenstates as an orthonormal basis
    `(es, ws)`, refusing to return one that cuts a degenerate multiplet.

    One extra state is computed so the gap *across* the cut can be seen.
    If `es[n] - es[n-1]` is negligible against the manifold's own energy
    spread, the retained states are not a symmetry-invariant subspace and
    any effective Hamiltonian fitted on them is meaningless -- so this
    raises ValueError, naming the nearest values of `n` that do cut
    cleanly, rather than returning numbers that cannot be trusted.
    """
    if n < 1: raise ValueError("get_manifold: n must be at least 1")
    try:
        es, ws = self.get_excited_states(n=n+1, mode=mode)
    except Exception: # cannot reach n+1 states (e.g. n+1 > Hilbert dim)
        import warnings
        warnings.warn("effectivehamiltonian: could not compute state n+1="
                      + str(n+1) + ", so the manifold at n=" + str(n)
                      + " could not be checked for a degeneracy cut -- the "
                      "returned couplings are only meaningful if that cut "
                      "happens to be clean")
        es, ws = self.get_excited_states(n=n, mode=mode)
        return np.array(es).real, list(ws)[:n]
    es = np.array(es).real
    ws = list(ws)
    if len(es) > n: # the extra state is available: check the cut
        scale = max(1.0, abs(es[n-1]-es[0]))
        gap = es[n]-es[n-1]
        if gap < degeneracy_tol*scale:
            safe = [k for k in range(1, len(es))
                    if es[k]-es[k-1] >= degeneracy_tol*scale]
            raise ValueError(
                "effectivehamiltonian: n=%d cuts a degenerate multiplet "
                "(E[%d]=%.12g and E[%d]=%.12g differ by %.3g). The retained "
                "states are then not a symmetry-invariant subspace, and the "
                "fitted couplings would be meaningless. Values of n that do "
                "cut cleanly here: %s"
                % (n, n-1, es[n-1], n, es[n], gap,
                   safe if safe else "none among the states computed"))
    return es[:n], orthonormalize(ws[:n])


def orthonormalize(ws, tol=ORTHONORMALITY_TOL):
    """Return an orthonormal basis spanning the same space as `ws`.

    `get_representation` assumes <psi_i|psi_j> = delta_ij. ED satisfies
    that to machine precision, but DMRG's excited states come from the
    overlap-penalty method and are only approximately orthogonal (~1e-9 on
    a small chain, and worse as system size or bond dimension grow), so
    the Gram matrix is checked and, when it needs to be, the basis is
    Lowdin-orthonormalized (S^-1/2, the orthonormal basis closest to the
    original one).

    The returned objects are linear combinations of the inputs, which
    every backend's MPS supports through __add__/__mul__.
    """
    ne = len(ws)
    S = np.zeros((ne, ne), dtype=np.complex128)
    for i in range(ne):
        for j in range(ne):
            S[i, j] = ws[i].overlap(ws[j])
    if np.max(np.abs(S-np.identity(ne))) < tol:
        return list(ws) # already orthonormal, nothing to do
    evals, evecs = np.linalg.eigh(S)
    if np.min(evals) <= 0.0:
        raise ValueError("effectivehamiltonian: the eigenstates are linearly "
                         "dependent (Gram matrix has a non-positive "
                         "eigenvalue %.3g) and cannot span a manifold"
                         % np.min(evals))
    sinvhalf = evecs@np.diag(1.0/np.sqrt(evals))@evecs.conj().T
    # built term by term starting from the first one, rather than with
    # sum(): the ED backend's State has no __radd__, so summing from the
    # implicit integer 0 raises TypeError
    out = []
    for i in range(ne):
        w = sinvhalf[0, i]*ws[0]
        for j in range(1, ne): w = w + sinvhalf[j, i]*ws[j]
        out.append(w)
    return out


def get_effective_hamiltonian_coefficients(self, mode="DMRG", n=4,
        tol=1e-4, operators=None):
    """Fit the projected Hamiltonian with a caller-supplied list of
    operators, and return one coefficient per operator (in the same
    order), zero where the fit put no weight.

    `tol` is the magnitude below which a fitted coefficient is reported as
    exactly zero.
    """
    if operators is None:
        raise ValueError("get_effective_hamiltonian_coefficients requires an "
                         "explicit list of operators to fit the projected "
                         "Hamiltonian with (operators=...)")
    es, ws = get_manifold(self, n, mode=mode)
    ops = {i: operators[i] for i in range(len(operators))}
    opm = {key: get_representation(ws, ops[key]) for key in ops}
    h = get_representation(ws, self.hamiltonian)
    opm["Id"] = np.identity(len(es)) # the overall energy offset
    coef = fit_matrix(h, opm, cutoff=tol)
    coef.pop("Id", None)
    # one entry per requested operator: a coefficient that fell below the
    # cutoff must come back as 0.0, not shorten the list and silently
    # misalign it against `operators`
    return [coef.get(i, 0.0) for i in range(len(operators))]


def get_effective_hamiltonian_couplings(self, mode="DMRG", n=4,
        method="single", tol=1e-4, operators=None):
    """Fit the projected Hamiltonian with the *pairwise products* of a set
    of single-site operators, and return a dict keyed by the pair.

    `operators`: dict of {name: MultiOperator} to build the products from;
    defaults to every Sx/Sy/Sz of the chain (`get_projection_operators`).

    `method` picks which operator the products actually mean, and the two
    are genuinely different operators rather than one being an
    optimization of the other:

    * `"single"` builds `P A P B P` -- the product of the two *projected*
      operators. Cheap (one projection per operator, then matrix
      products), and equal to the below only when the manifold happens to
      be closed under `A`.
    * `"full"` builds `P (A B) P` -- the projection of the product. One
      projection per *pair*, so quadratically more expensive, but it is
      the operator the notation suggests.

    Note that a product of two Hermitian matrices is not itself Hermitian,
    so with either method the fit basis is not Hermitian even though the
    target is; the anti-Hermitian parts have to cancel among themselves.

    `tol` is the magnitude below which a fitted coefficient is dropped.
    """
    es, ws = get_manifold(self, n, mode=mode)
    op = get_projection_operators(self) if operators is None else operators
    opm = dict()
    if method == "single":
        for key in op: opm[key] = get_representation(ws, op[key])
    elif method != "full":
        raise ValueError("get_effective_hamiltonian_couplings: method must be "
                         "'single' or 'full', got "+repr(method))
    ops = {"Id": np.identity(len(es))}
    for key1 in op: # use just the quadratic operators
        for key2 in op:
            if method == "single":
                m = opm[key1]@opm[key2]
            else:
                m = get_representation(ws, op[key1]*op[key2])
            if acceptable_matrix(m, ops): # if the matrix can be accepted
                ops[(key1, key2)] = m
    h = get_representation(ws, self.hamiltonian)
    coef = fit_matrix(h, ops, cutoff=tol)
    coef.pop("Id", None)
    return coef


def get_effective_hamiltonian(self, tol=1e-3, **kwargs):
    """Compute the effective Hamiltonian and return its latex form"""
    coef = get_effective_hamiltonian_couplings(self, **kwargs)
    return dict2latex(coef, tol=tol) # return the latex form


def dict2latex(d, tol=1e-3):
    """Transform the dictionary of couplings into a latex display block,
    with every coupling normalized by the largest one."""
    if len(d) < 1: return r"H = 0" + "\n"
    cs = [d[key] for key in d] # coefficients
    cmax = np.round(np.max(np.abs(cs)), 3) # maximum value
    if cmax <= 0.0: # every coupling rounds to zero: nothing to normalize by
        return r"H = 0" + "\n"
    terms = [] # the surviving terms, with their own signs
    for key in d: # loop
        c = np.round(d[key]/cmax, 3) # round the number
        if np.abs(c) < tol: continue
        names = key if isinstance(key, tuple) else (str(key),)
        terms.append((c, "  ".join(names)))
    if len(terms) < 1: return r"H = 0" + "\n"
    # the sign becomes the separator, so a negative coupling reads
    # "- 0.587 ..." rather than "+ -0.587 ...", and the last term no
    # longer trails a dangling "+"; the prefactor lives inside the math
    # block rather than in front of it
    body = ("-  " if terms[0][0] < 0 else "") + str(abs(terms[0][0])) \
        + "  " + terms[0][1]
    for c, name in terms[1:]:
        body += "\n" + ("-" if c < 0 else "+") + "  " + str(abs(c)) \
            + "  " + name
    return r"\[" + "\n" + "H = " + str(cmax) + r" \left(" + "\n" \
        + body + "\n" + r"\right)" + "\n" + r"\]" + "\n"


def get_projection_operators(self):
    """The single-site spin operators of every site of the chain, keyed by
    a latex-ish name.

    The site count comes from the length of the chain's own Sx list, not
    from self.ns: for Spinful_Fermionic_Chain those differ (self.ns counts
    the two interleaved spinless sites per orbital, while Sx has one entry
    per orbital), and hardcoding self.ns//2 -- as this used to -- silently
    covered only the first half of the sites of every *other* class,
    including Spin_Chain.
    """
    op = dict() # dictionary
    for i in range(len(self.Sx)):
        op["S^x_"+str(i+1)] = self.Sx[i]
        op["S^y_"+str(i+1)] = self.Sy[i]
        op["S^z_"+str(i+1)] = self.Sz[i]
    return op # return operators


def get_representation(ws, op):
    """Return the representation of a certain operator"""
    ne = len(ws)
    h = np.zeros((ne, ne), dtype=np.complex128)
    for i in range(ne):
        for j in range(ne):
            h[i, j] = ws[i].overlap(op*ws[j])
    return h # return matrix


def acceptable_matrix(m, ops):
    """Check if it is ok to keep this matrix.

    A cheap pre-filter only: it drops matrices that are zero or nearly
    parallel to one already kept, which keeps the fit's design matrix
    small. It is *not* what makes the fit well posed -- pairwise
    near-orthogonality does not imply linear independence -- that is
    fit_matrix's own rank check.
    """
    if np.sum(np.abs(m)) < 1e-7: return False
    v = matrix2vector(m)
    for key in ops: # loop over the other matrices
        o = ops[key] # get the matrix
        vo = matrix2vector(o) # convert to vector
        proj = v.dot(vo)/(np.sqrt(v.dot(v))*np.sqrt(vo.dot(vo)))
        if np.abs(proj) > 0.98:
            return False
    return True


def matrix2vector(m):
    """Flatten a complex matrix into a real vector, stacking real and
    imaginary parts. The dot product of two such vectors is the real
    Frobenius inner product Re Tr(A^dag B)."""
    n = m.shape[0] # get the dimension
    v = np.zeros(2*n**2) # to a vector
    v[0:n**2] = m.reshape(n**2).real
    v[n**2:2*n**2] = m.reshape(n**2).imag
    return v


def fit_matrix(h, d, cutoff=1e-4):
    """Fit a matrix with a dictionary of matrices: find real x minimizing
    ||h - sum_i x_i d_i||, and return {key: x_i} for the coefficients
    above `cutoff`.

    This is an ordinary linear least-squares problem, and is solved as
    one. It used to be handed to scipy.optimize.minimize from an unseeded
    np.random.random start, which cost about eight digits of accuracy
    (measured residual 3.5e-6 against lstsq's 1.4e-14 on the same system)
    and made the result irreproducible: two successive calls differed by
    ~2.6e-5, enough for a coupling near `cutoff` to be returned by one
    call and dropped by the next, changing the key set of the result.

    Raises ValueError if the candidate matrices are linearly dependent on
    this manifold. That case has no unique answer -- the coefficients can
    be shifted by anything in the null space without changing the fit --
    so returning one arbitrary representative would be reporting a
    specific set of couplings that the data does not actually determine.
    """
    keys = list(d) # fix an order; dicts preserve insertion order
    D = np.array([matrix2vector(d[key]) for key in keys]).T
    y = matrix2vector(h)
    x, _, rank, sv = np.linalg.lstsq(D, y, rcond=None)
    if rank < len(keys):
        raise ValueError(
            "effectivehamiltonian: the %d candidate operators are linearly "
            "dependent on this manifold (rank %d, so %d null directions), so "
            "no unique set of couplings exists -- the fit could be shifted by "
            "any null-space vector without changing its residual. Use fewer "
            "or more independent operators, or a larger manifold (bigger n)."
            % (len(keys), rank, len(keys)-rank))
    out = dict()
    for i, key in enumerate(keys): # loop over the operators
        if np.abs(x[i]) > cutoff:
            out[key] = x[i]
    return out # return the coefficients
