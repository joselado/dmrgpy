"""Non-Hermitian DMRG on the live Julia backend.

Both entry points here return one *attempt* -- a single
(energy,psil,psir) triple from a fresh random start. The retry loop and
the two-sided eigen-residual certificate that decides whether an attempt
converged live one level up, in the backend-agnostic nhdmrg.py
(nhdmrg()/nhdmrg_generalized()): those are built purely on generic
MultiOperator*MPS algebra, which this backend's MPS already supports, so
they are shared rather than reimplemented here.

Unlike the other three backends, the local non-Hermitian sweep itself is
not dmrgpy's own port of ITensorNHDMRG.jl -- it *is* ITensorNHDMRG.jl
(declared in ../juliapkg.json, already used by get_gs.jl's
get_gs_nhdmrg for plain non-Hermitian ground states). So the two
algorithm knobs the ported backends expose, krylovdim and restarts, do
not map one-to-one: see nhdmrg_attempt's own note.
"""

from .generalized import metric_guard
from .juliasession import Main as Mainjl
from .mpo import MPO
from .mps import MPS, random_mps
from .groundstate import NH_alg, NH_biorthoalg


def julia_random_mps(self):
    """A fresh random MPS living in the Julia session.

    Deliberately mpsjulialive.mps.random_mps rather than the chain-level
    Many_Body_Chain.random_state(): that one honors self.mode, so on a
    chain with mode="ED" it returns an edtk State, which has no .jlmps
    and blows up with an opaque AttributeError several frames later
    inside the juliacall call below. The session backends never hit this
    because they build their random start inside their own C++/Python
    session rather than through the chain-level helper -- so routing
    through random_state() here would have made julia_live the only
    backend where mode="ED" turns nhdmrg() into an AttributeError instead
    of just running the DMRG solver, which is what every backend does
    (NH-DMRG has no ED implementation for mode= to select).
    """
    return random_mps(self)


def nhdmrg_attempt(self, H, krylovdim=20, restarts=2):
    """One NH-DMRG run from a fresh random start. Returns
    (energy,psil,psir) in dmrgpy's own biorthogonal convention
    (H^dagger|psil>=conj(energy)|psil>, <psil|psir>=1) -- which is *not*
    ITensorNHDMRG's own, see nhdmrg.jl's header for the conjugation that
    get_nhdmrg_pair applies.

    krylovdim maps onto ITensorNHDMRG's eigsolve_krylovdim, but restarts
    does *not* map onto its eigsolve_maxiter: the two count different
    things (dmrgpy's "restarts" is how many times the ported backends
    restart their own Arnoldi from the current best Ritz vector, whereas
    eigsolve_maxiter bounds KrylovKit's internal restarts of a solve that
    is already differently structured). Rather than pretend they're
    interchangeable, restarts is accepted and ignored here, and
    eigsolve_maxiter keeps get_gs.jl's own default -- the outer DMRG
    sweeps (self.nsweeps) are what actually converge the run on every
    backend anyway, which is why nhdmrg.py keeps both knobs small.
    """
    Hop = MPO(H, MBO=self)
    psi0 = julia_random_mps(self)
    with metric_guard(): # nh_biorthogonal_pair's own collapse guard
        e, psil, psir = Mainjl.get_nhdmrg_pair(Hop.jlmpo, psi0.jlmps,
                nsweeps=self.nsweeps, cutoff=self.cutoff, maxm=self.maxm,
                alg=NH_alg, biorthoalg=NH_biorthoalg,
                eigsolve_krylovdim=int(krylovdim))
    return complex(e), MPS(psil, MBO=self), MPS(psir, MBO=self)


def nhdmrg_generalized_attempt(self, H, A, krylovdim=20, restarts=2,
        lam0=None):
    """One non-Hermitian generalized-eigenproblem run
    (H|psi_R>=lambda*A|psi_R>) from a fresh random start, via
    generalized.jl's get_gs_generalized_nhdmrg. Returns
    (lam,psil,psir), same conventions as nhdmrg_attempt above; restarts
    is accepted and ignored for the same reason.

    lam0 (optional starting lambda estimate) is passed through as a
    complex with a NaN real part when unset -- generalized.jl's own
    sentinel, matching Chain::nhdmrg_generalized's pybind11 binding.
    """
    Hop = MPO(H, MBO=self)
    Aop = MPO(A, MBO=self)
    psi0 = julia_random_mps(self)
    jl_lam0 = complex(float("nan"), 0.0) if lam0 is None else complex(lam0)
    with metric_guard(): # near-null-space guard -> RuntimeError, which
                          # nhdmrg.py's retry loop treats as a bad draw
        lam, wfl, wfr = Mainjl.get_gs_generalized_nhdmrg(Hop.jlmpo,
                Aop.jlmpo, psi0.jlmps, nsweeps=self.nsweeps,
                cutoff=self.cutoff, maxm=self.maxm,
                lam0=jl_lam0, alg=NH_alg, biorthoalg=NH_biorthoalg,
                eigsolve_krylovdim=int(krylovdim))
    return complex(lam), MPS(wfl, MBO=self), MPS(wfr, MBO=self)
