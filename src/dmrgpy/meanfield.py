# routines to perform mean field calculations
import numpy as np

from . import multioperator

# Component order used throughout this module: the index a in S^a is 0/1/2
# for Sx/Sy/Sz, matching Spin_Chain.get_magnetization()'s own row order.
_COMPONENTS = ("Sx","Sy","Sz")
_COMPONENT_INDEX = {n:k for k,n in enumerate(_COMPONENTS)}

# Coefficient magnitude below which a decomposed term is treated as absent.
# Same spirit as multioperator.clean_threshold, and applied for the same
# reason: a Hamiltonian written as a sum of SS(i,j) products carries exact
# zeros for the components that cancel, and rebuilding those as operator
# terms would triple the size of the mean-field Hamiltonian for nothing.
_tol = 1e-10


def decompose_spin_hamiltonian(sc):
    """Split a spin chain's Hamiltonian into its one-site and two-site
    parts, returned as plain arrays:

      b[i,a]      coefficient of S^a_i          (the model's own fields)
      J[i,j,a,b]  coefficient of S^a_i S^b_j    (the exchange to decouple)
      const       coefficient of the identity   (a constant energy offset)

    This reads `sc.hamiltonian`, the MultiOperator the chain was actually
    given by `set_hamiltonian()`. It replaces the old `sc.exchange` list,
    which `Spin_Chain.set_exchange()` used to populate and which has been
    the integer 0 -- and this whole module therefore broken -- ever since
    that builder was removed in favour of writing Hamiltonians out with
    `SS(i,j)`.

    Note J is indexed by the *ordered* pair: a term written S^x_0 S^x_1
    lands in J[0,1,0,0] and one written S^x_1 S^x_0 in J[1,0,0,0]. The
    Weiss field below sums both orders, so either spelling (or a mixture)
    gives the same answer.

    Identity factors are dropped before a term is classified: building a
    Hamiltonian as `h = 0` followed by `h = h + ...` leaves a literal
    `[0.0, ['Id', 1]]` term behind, so "Id" is an ordinary thing to meet
    here and a term made only of identities is just a constant.

    Raises ValueError for anything a mean-field decoupling of a spin
    Hamiltonian is not defined for: a non-spin operator, a product of
    three or more sites, or two factors on the same site (S^a_i S^b_i is a
    one-site operator, not an exchange bond -- decoupling it would replace
    an exact on-site term with an approximation for no reason).
    """
    h = sc.hamiltonian
    if h is None:
        raise ValueError("meanfield: this chain has no Hamiltonian yet; "
                         "call set_hamiltonian() first")
    ns = sc.ns
    b = np.zeros((ns,3),dtype=np.complex128)
    J = np.zeros((ns,ns,3,3),dtype=np.complex128)
    const = 0.0
    for term in h.op: # [coefficient, [name, site], [name, site], ...]
        c,factors = term[0],term[1:]
        factors = [f for f in factors if f[0]!="Id"] # identities carry nothing
        for name,i in factors:
            if name not in _COMPONENT_INDEX:
                raise ValueError(
                    "meanfield: operator %r is not a spin component; only "
                    "Hamiltonians built from Sx/Sy/Sz can be mean-field "
                    "decoupled by this routine"%name)
        if len(factors)==0:
            const += c
        elif len(factors)==1:
            (name,i), = factors
            b[i,_COMPONENT_INDEX[name]] += c
        elif len(factors)==2:
            (n0,i0),(n1,i1) = factors
            if i0==i1:
                raise ValueError(
                    "meanfield: term %s_%d %s_%d acts twice on site %d; "
                    "that is a one-site operator, not an exchange bond, "
                    "and has no mean-field decoupling"%(n0,i0,n1,i1,i0))
            J[i0,i1,_COMPONENT_INDEX[n0],_COMPONENT_INDEX[n1]] += c
        else:
            raise ValueError(
                "meanfield: term acting on %d sites; only one- and "
                "two-site terms can be mean-field decoupled"%len(factors))
    return b,J,const


def _weiss_field(J,m):
    """The Weiss field h[i,a] = sum_{j,b} (J[i,j,a,b] + J[j,i,b,a]) m[j,b].

    This is the standard decoupling S^a_i S^b_j -> S^a_i <S^b_j> +
    <S^a_i> S^b_j (minus the constant), summed over both orders so that
    each site picks up the field from its partner whichever way round the
    bond was written -- the general-coupling version of the old code's
    symmetric `fnew[c.i] += g.m[c.j]; fnew[c.j] += g.m[c.i]` pair.
    """
    return np.einsum("ijab,jb->ia",J,m) + np.einsum("jiba,jb->ia",J,m)


def _mean_field_hamiltonian(sc,bonds,b,const,hmf,p):
    """H_MF = p*H_exchange + H_onsite + (1-p)*sum_i h_i.S_i + const.

    The model's own one-site terms are kept in full: they are an external
    field, not something being decoupled. At p=1 the Weiss field drops out
    and this is exactly the Hamiltonian that came in; at p=0 the exchange
    drops out and this is pure mean-field theory.
    """
    Si = [sc.Sx,sc.Sy,sc.Sz]
    h = 0
    if abs(p)>_tol: # many-body exchange, kept at strength p
        for (i,j,a,bb,g) in bonds: h = h + p*g*Si[a][i]*Si[bb][j]
    ns = sc.ns
    for i in range(ns): # the model's own fields, plus the Weiss field
        for a in range(3):
            c = b[i,a] + (1.-p)*hmf[i,a]
            if abs(c)>_tol: h = h + c*Si[a][i]
    if abs(const)>_tol: h = h + const*sc.Id # constant energy offset
    if not isinstance(h,multioperator.MultiOperator):
        # h is still the integer 0: every coefficient fell below _tol, so
        # there is nothing to solve and set_hamiltonian would be handed a
        # bare int.
        raise ValueError("meanfield: the mean-field Hamiltonian is empty "
                         "(every coefficient vanished)")
    return h


def spinchain_meanfield(sc,p=0.0,mix=0.9,m0=None,maxerror=1e-06,
        maxite=1000,mixmode="default",**kwargs):
    """Mean field calculation of a spin chain,
    p controls the degree of many body correlations.

    The chain's Hamiltonian is read as it was written with
    `set_hamiltonian()` (see `decompose_spin_hamiltonian`), so any
    Hamiltonian built out of Sx/Sy/Sz one- and two-site terms works --
    including one whose exchange is anisotropic or whose bonds are not
    nearest-neighbour. Extra keyword arguments are forwarded to
    `gs_energy`/`get_magnetization`, so `mode="ED"` picks the ED solver.
    (The mixing scheme is `mixmode=`, not `mode=`: `mode` means the
    DMRG/ED solver everywhere else in this library, and this routine
    forwards it there.)

    Returns the converged chain, whose Hamiltonian is the last mean-field
    one and whose ground state is its solution.
    """
    sc0 = sc.copy() # copy the spin chain object
    b,J,const = decompose_spin_hamiltonian(sc0) # read the Hamiltonian once
    # The exchange never changes across iterations, so flatten it to the
    # nonzero entries once instead of walking a mostly-empty ns^2 x 9
    # array per iteration.
    bonds = [(i,j,a,bb,J[i,j,a,bb])
             for i in range(sc0.ns) for j in range(sc0.ns)
             for a in range(3) for bb in range(3)
             if abs(J[i,j,a,bb])>_tol]
    if m0 is None:
        mold = np.array([np.random.random(3) for i in range(sc0.ns)])
    else: mold = np.array(m0)
    def get_new(mold):
        """Return new magnetization"""
        sc = sc0.copy() # copy the initial Hamiltonian
        hmf = _weiss_field(J,np.array(mold,dtype=np.complex128))
        sc.set_hamiltonian(_mean_field_hamiltonian(sc,bonds,b,const,hmf,p))
        sc.gs_energy(**kwargs) # perform calculation
        mnew = np.array(sc.get_magnetization(**kwargs)).transpose()
        return mnew,sc # return new magnetization
    if mixmode=="default": # mixing scheme
      sc = sc0
      for ite in range(maxite):
          mnew,sc = get_new(mold) # new magnetization
          error = np.sum(np.abs(mnew-mold)) # error, before mixing
          mold = mix*mnew + (1.0-mix)*mold # mix for the next iteration
          if error<maxerror: break
          print("Error",error)
      else:
          print("WARNING: mean field did not converge in",maxite,
                "iterations, last error",error)
    elif mixmode=="broyden": # mixing scheme
        from scipy.optimize import broyden1
        def fopt(mold): # function to optimize
          mnew,sc = get_new(mold) # new magnetization
          return mnew - mold
        try: mnew = broyden1(fopt,mold,maxiter=20)
        except: return spinchain_meanfield(sc0,p=p,mix=mix,m0=mold,
                maxerror=maxerror,maxite=maxite,**kwargs)
        print("Done")
        mnew,sc = get_new(mnew) # new magnetization
    else:
        raise ValueError("meanfield: unrecognized mixmode %r; expected "
                         "'default' or 'broyden'"%mixmode)

    print("Done")
    return sc # return the spin chain object
