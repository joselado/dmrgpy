import numpy as np
from . import mps


def get_excited_states_dmrg(self,n=2,noise=0.0,scale=10.0):
    """
    Excited states via the in-process pybind11 extension
    (mpscpp2/chain_session.h's Chain::excited_states). The "noise" kwarg is
    intentionally not passed through: in the old file-based backend it was
    written to tasks.in under the key "noise", but get_sweeps.h actually
    read the (differently spelled) key "moise" for the Sweeps object, so
    that kwarg never actually affected anything even before this backend
    was removed -- not reproducing a pre-existing no-op.
    """
    if self.itensor_version=="julia_live":
        from .mpsjulialive import excited as exjl
        return exjl.get_excited_states_dmrg(self,n=n,noise=noise,scale=scale)
    wf0 = self.get_gs()
    self._session.set_sweep_params(self.maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_verbose(self.verbose)
    self._session.set_mpomaxm(max(self.maxm,self.mpomaxm))
    energies,fluctuations,handles = self._session.excited_states(
            n,scale,self.excited_gram_schmidt)
    wfs = [mps.MPS(MBO=self,cpp_handle=h).copy() for h in handles]
    return np.array(energies),wfs


def get_excited(*args,**kwargs):
    """Return the excited state energies"""
    (es,ws) = get_excited_states(*args,**kwargs)
    return np.array(es)



def remove_none(ws):
    wout = []
    for w in ws:
        if w is not None: wout.append(w)
    return wout



def get_excited_states(self,n=2,purify=True,**kwargs):
    """Excited states"""
    H = self.hamiltonian # make a check
    if not self.is_hermitian(H): # non-Hermitian Hamiltonians
        if n==1 and self.itensor_version in (2,3,"python"):
            # ground state only: use the same NH-DMRG solve gs_energy()
            # routes to (groundstate.py's non-Hermitian branch), so the
            # two entry points cannot return different eigenpair members
            # of a degenerate multiplet or differently-normalized states
            e0 = self.gs_energy(**kwargs)
            return (np.array([e0]),[self.wf0.copy()])
#        print("Non-Hermitian Hamiltonian, using Arnoldi method")
        return excited_states_non_hermitian(self,n=n,**kwargs)
    if n==1: # workaround for just the ground state
        w = self.get_gs(**kwargs)
        e0 = self.gs_energy(**kwargs)
        return ([e0],[w]) # return
    if not purify: # just compute excited states
        return get_excited_states_dmrg(self,n=n,**kwargs) # compute 
    else: # purify the states
        es,ws = get_excited_states_dmrg(self,n=n+2,**kwargs) # compute
        ws = gram_smith(ws) # orthogonalize the MPS
        ws = remove_none(ws) # remove None wavefunctions
        ne = len(ws) # number of states
        # rediagonalize the Hamiltonian in this subspace once, taking
        # both refined eigenvalues and eigenvectors from the same
        # O(ne^2) representation matrix (this used to be built twice:
        # once here just for eigvalsh's eigenvalues, then again from
        # scratch inside rediagonalize for the eigenvectors)
        from .algebra.arnolditk import rediagonalize
        es,ws = rediagonalize(self.hamiltonian,ws) # rediagonalize
        if ne<n: n = ne # redefine
        return (es[0:n],ws[0:n])


from .algebra.arnolditk import gram_smith


def excited_states_non_hermitian(self,n=3,recursive=True,
        maxit=40,nkry_min =3,nkry_max =10,
     **kwargs):
    # recursive=True (one state at a time, each deflated against the
    # previous ones) is the default rather than a single simultaneous
    # (n>1) Arnoldi search: a shared Krylov subspace across several
    # requested eigenvalues can fail to resolve a (near-)degenerate pair
    # -- confirmed on a non-Hermitian chain with a 2-fold-degenerate
    # excited state, where the simultaneous search converged to two
    # spurious values instead. The one-at-a-time search is slower but
    # reliable, since each state gets its own warm start explicitly
    # deflated against the others.
    from . import mpsalgebra
    (es,wf) = mpsalgebra.mpsarnoldi(self,self.hamiltonian,mode="GS",
                nwf=n, # number of wavefunctions
                recursive_arnoldi=recursive, # recursive?
                maxit = maxit, # max number of Arnoldi iterations
                nkry_min = nkry_min, # minimum number of Krylov vectors
                nkry_max = nkry_max, # maximum number of Krylov vectors
                **kwargs)
    from .algebra.algebra import sorteigen
    es,wf = sorteigen(es,wf)
    return (es,wf)


