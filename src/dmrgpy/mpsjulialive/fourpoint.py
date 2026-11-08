"""Four-point correlator tensor <Cdag_i C_j Cdag_k C_l> on julia_live,
computed in Julia by mpsjulialive/fourpoint.jl's four_point_sweep (one
walk over the sites that shares left environments between the tuples
with a common prefix, Jordan-Wigner strings written site by site)."""
import numpy as np


def four_correlation_tensor_sweep(wf, accelerate=True):
    """ctmode="sweep" for an itensor_version="julia_live" wavefunction.

    Plain spinless Fermion sites only (Fermionic_Chain and the interleaved
    Spinful_Fermionic_Chain), the same scope as the v3 and "python" sweeps.
    accelerate is accepted for the common signature: every entry is a leaf
    of the same walk, so there is no per-tuple work to skip. Returns the raw
    <wf|O|wf>, without normalizing wf."""
    from ..entropytk.correlationentropy import _julia_sweep_applies
    from .juliasession import Main as Mainjl, to_julia_strvec
    if not _julia_sweep_applies(wf.MBO):
        raise NotImplementedError("ctmode=\"sweep\" on julia_live needs a "
                "chain of plain spinless fermion sites with one C/Cdag "
                "per site; use ctmode=\"explicit\" instead")
    names = to_julia_strvec(["Cdag", "C", "Cdag", "C"])
    out = Mainjl.four_point_sweep(wf.jlmps, names)
    return np.array(out, dtype=np.complex128)
