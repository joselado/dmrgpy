import numpy as np

from .. import multioperator


def vev(MBO,MO,wf=None,npow=1,T=0.):
    """<wf|MO^npow|wf> on julia_live, the counterpart of vev.py's
    multi_vev on the session backends.

    The state is the chain's ground state when wf is None, and a unit-norm
    copy of wf otherwise (groundstate.unit_copy), so a state given as 2*s
    reads as s. For npow>1 the MPO is applied npow-1 times, each
    application truncated to the chain's maxm and cutoff, as on the session
    backends; the last factor is contracted exactly, <a|M|b> with no
    truncation (gs_extra.jl's gs_extra_aMb). No **kwargs: an unknown
    keyword raises TypeError."""
    if T!=0.:
        raise NotImplementedError("vev with T>0 has no DMRG route: use "
                "mode=\"ED\" or metts_vev()")
    MO = multioperator.obj2MO(MO,name="vev_multioperator")
    if npow==0: return 1.0
    npow = int(npow)
    if npow<0: raise ValueError("vev: npow must be a non-negative integer")
    from .. import groundstate
    if wf is None: wf = MBO.get_gs() # normalized already
    else: wf = groundstate.unit_copy(wf)
    from .mpo import MPO
    from .juliasession import Main as Mainjl
    M = MPO(MO,MBO=MBO)
    phi = wf.jlmps
    for i in range(npow-1):
        phi = Mainjl.gs_extra_apply(M.jlmpo,phi,int(MBO.maxm),
                float(MBO.cutoff))
    return complex(Mainjl.gs_extra_aMb(wf.jlmps,M.jlmpo,phi))


# the same cap as vev.py's: the result's bond dimension is bounded by the
# MPO's times the state's, far below this, so only the cutoff acts
_UNCAPPED = 1000000


def energy_variance(MBO,h,wf=None):
    """<wf|(H-<H>)^2|wf> on julia_live, as ||(H-<H>)|wf>||^2: the MPO of
    h-<h> applied once with the maxm cap lifted, so only the relative
    cutoff acts, and the result contracted exactly against the state.
    vev.py's energy_variance is the rule; this one has no session to
    reset afterwards, only the chain's maxm."""
    if wf is None: wf = MBO.get_gs()
    e = vev(MBO,h,wf=wf)
    dh = h - e*multioperator.identity()
    old = MBO.maxm
    MBO.maxm = _UNCAPPED
    try:
        out = vev(MBO,dh,wf=wf,npow=2)
    finally:
        MBO.maxm = old
    return out
