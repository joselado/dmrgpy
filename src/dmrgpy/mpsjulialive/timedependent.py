import numpy as np

from .. import operatornames
from .. import multioperator
from .. timedependent import _is_tebd_nn_error
from .mpo import MPO, text_mpo


def _tebd_or_tdvp(tebd_call,tdvp_call):
    """Julia-backend counterpart of timedependent.py's own
    _tebd_or_tdvp(): backs tevol_method="AUTO" by trying `tebd_call`
    first and transparently retrying as `tdvp_call` if tebd.jl's
    bond_hamiltonians() rejects self.hamiltonian for not being strictly
    nearest-neighbor. juliacall surfaces every Julia-side `error(...)` as
    a JuliaError regardless of which Julia function raised it, so this
    reuses the same message-text check (_is_tebd_nn_error) as the
    C++/pure-Python backends to avoid swallowing an unrelated Julia
    failure (a real ITensorMPS/KrylovKit bug, say) as if it meant "not
    nearest-neighbor"."""
    from juliacall import JuliaError
    try:
        return tebd_call()
    except JuliaError as exc:
        if not _is_tebd_nn_error(exc): raise
        return tdvp_call()


def _check_tevol_method(self,methods=("TDVP","TDVP_GSE")):
    """This backend implements self.tevol_method in `methods` only for
    whichever operation is calling this (default: ("TDVP","TDVP_GSE"), the
    set advance_complex_time_step supports; evolution_dmrg_DC/
    evolve_and_measure_dmrg below pass a wider set that also includes
    "TEBD" -- TDZ/complex-time evolution has no TEBD counterpart on any
    backend, see tdz.py's own _advance_complex_time_step, which never
    branches on "TEBD" either). Raised rather than silently falling back to
    plain TDVP: a caller who set tevol_method="TEBD" asked for a specific
    integrator, and quietly running a different one would make a
    backend-comparison script silently compare the wrong things."""
    if self.tevol_method not in methods:
        raise NotImplementedError(
            "itensor_version='julia_live' implements tevol_method "
            +", ".join(repr(m) for m in methods)+" for this operation "
            "(got "+repr(self.tevol_method)+")")


def evolution_dmrg_DC(self,name="XX",nt=10000,dt=0.1,**kwargs):
    """Real-time quench dynamical correlator on the Julia backend, via
    native TDVP (mpsjulialive/tdvp.jl's quench_tdvp -- the whole nt-step
    trajectory runs in one Julia call, same design as kpm.jl), its one-site
    + global-subspace-expansion variant (quench_tdvp_gse) when
    self.tevol_method=="TDVP_GSE", or 2nd-order-Trotter TEBD
    (mpsjulialive/tebd.jl's quench_tebd) when self.tevol_method=="TEBD"
    (self.tevol_method=="AUTO" tries this first and falls back to plain
    "TDVP" -- via _tebd_or_tdvp() above -- if tebd.jl rejects
    self.hamiltonian for not being nearest-neighbor). Mirrors
    timedependent.py::evolution_dmrg_DC's TDVP/TDVP_GSE/TEBD/AUTO
    branches (itensor_version=3/"python"); the legacy MPO-Taylor path has
    no julia_live implementation at all."""
    _check_tevol_method(self,methods=("TDVP","TDVP_GSE","TEBD","AUTO"))
    name = operatornames.str2MO(self,name,**kwargs)
    name[0] = name[0].get_dagger()
    A,B = name[0],name[1]
    wf0 = self.get_gs() # also sets self.e0
    H = self.hamiltonian
    # shift H by the ground-state energy so psi1's trivial overall phase
    # exp(-i*e0*t) cancels in the <psi2|psi1(t)> correlator, mirroring
    # pyitensor.chain.Chain.quench_tdvp()/quench()
    Hshift_MO = H-self.e0*multioperator.identity()
    A1 = MPO(A,MBO=self)
    A2 = MPO(B,MBO=self)
    from .juliasession import Main as Mainjl, to_julia_strvec
    if self.tevol_method=="TEBD":
        # text_mpo() (mpo.py, shared with the MPO path) serializes the
        # *raw*, non-Jordan-Wigner-dressed term list ("C"/"Cdag" names) --
        # exactly what tebd.jl's read_terms()/resolve_term() need, and the
        # only form real ITensors.jl's builtin "Fermion" site type
        # understands natively (see tebd.jl's own header comment).
        hlines = to_julia_strvec(text_mpo(Hshift_MO))
        correlator,_wf = Mainjl.quench_tebd(hlines,self.jlsites,
                A1.jlmpo,A2.jlmpo,wf0.jlmps,int(nt),dt,self.cutoff,self.maxm)
    elif self.tevol_method=="AUTO":
        hlines = to_julia_strvec(text_mpo(Hshift_MO))
        correlator,_wf = _tebd_or_tdvp(
                lambda: Mainjl.quench_tebd(hlines,self.jlsites,
                    A1.jlmpo,A2.jlmpo,wf0.jlmps,int(nt),dt,self.cutoff,self.maxm),
                lambda: Mainjl.quench_tdvp(MPO(Hshift_MO,MBO=self).jlmpo,
                    A1.jlmpo,A2.jlmpo,wf0.jlmps,int(nt),dt,self.cutoff,self.maxm))
    else:
        Hshift = MPO(Hshift_MO,MBO=self)
        if self.tevol_method=="TDVP_GSE":
            correlator,_wf = Mainjl.quench_tdvp_gse(Hshift.jlmpo,
                    A1.jlmpo,A2.jlmpo,wf0.jlmps,int(nt),dt,self.cutoff,self.maxm,
                    int(self.tdvp_gse_sweeps),int(self.tdvp_gse_krylov_order),
                    self.tdvp_gse_cutoff)
        else:
            correlator,_wf = Mainjl.quench_tdvp(Hshift.jlmpo,
                    A1.jlmpo,A2.jlmpo,wf0.jlmps,int(nt),dt,self.cutoff,self.maxm)
    cs = np.array(correlator)
    ts = np.array([dt*ii for ii in range(int(nt))])
    return ts,cs.real-1j*cs.imag


def evolve_and_measure_dmrg(self,operator=None,nt=1000,h=None,
        dt=1e-2,wf=None,return_wf=False,**kwargs):
    """Real-time evolution + measurement on the Julia backend, via native
    TDVP (mpsjulialive/tdvp.jl's evolve_and_measure_tdvp), its one-site
    + global-subspace-expansion variant (evolve_and_measure_tdvp_gse) when
    self.tevol_method=="TDVP_GSE", or 2nd-order-Trotter TEBD
    (mpsjulialive/tebd.jl's evolve_and_measure_tebd) when
    self.tevol_method=="TEBD" (self.tevol_method=="AUTO" tries this first
    and falls back to plain "TDVP" -- via _tebd_or_tdvp() above -- if
    tebd.jl rejects `h` for not being nearest-neighbor). Mirrors
    timedependent.py::evolve_and_measure_dmrg's TDVP/TDVP_GSE/TEBD/AUTO
    branches."""
    _check_tevol_method(self,methods=("TDVP","TDVP_GSE","TEBD","AUTO"))
    if h is None: h = self.hamiltonian # Hamiltonian
    if wf is None: wf = self.wf0 # get ground state
    Aop = MPO(operator,MBO=self)
    from .juliasession import Main as Mainjl, to_julia_strvec
    if self.tevol_method=="TEBD":
        hlines = to_julia_strvec(text_mpo(h))
        correlator,wf_final_jl = Mainjl.evolve_and_measure_tebd(hlines,
                self.jlsites,Aop.jlmpo,wf.jlmps,int(nt),dt,self.cutoff,
                self.maxm)
    elif self.tevol_method=="AUTO":
        hlines = to_julia_strvec(text_mpo(h))
        correlator,wf_final_jl = _tebd_or_tdvp(
                lambda: Mainjl.evolve_and_measure_tebd(hlines,
                    self.jlsites,Aop.jlmpo,wf.jlmps,int(nt),dt,self.cutoff,
                    self.maxm),
                lambda: Mainjl.evolve_and_measure_tdvp(MPO(h,MBO=self).jlmpo,
                    Aop.jlmpo,wf.jlmps,int(nt),dt,self.cutoff,self.maxm))
    else:
        Hop = MPO(h,MBO=self)
        if self.tevol_method=="TDVP_GSE":
            correlator,wf_final_jl = Mainjl.evolve_and_measure_tdvp_gse(
                    Hop.jlmpo,Aop.jlmpo,wf.jlmps,int(nt),dt,self.cutoff,
                    self.maxm,int(self.tdvp_gse_sweeps),
                    int(self.tdvp_gse_krylov_order),self.tdvp_gse_cutoff)
        else:
            correlator,wf_final_jl = Mainjl.evolve_and_measure_tdvp(Hop.jlmpo,
                    Aop.jlmpo,wf.jlmps,int(nt),dt,self.cutoff,self.maxm)
    cs = np.array(correlator)
    ts = np.array([dt*ii for ii in range(int(nt))])
    if return_wf:
        from .mps import MPS
        wf_final = MPS(wf_final_jl,MBO=self)
        return ts,cs.real-1j*cs.imag,wf_final
    return ts,cs.real-1j*cs.imag


def advance_complex_time_step(self,Hop,wf,dz,do_gse=False):
    """One complex-time TDVP step exp(-i*dz*Hop) on the Julia backend --
    used by tdz.py's submode="TDZ" (arXiv:2311.10909). Reuses
    tdvp.jl's tdvp_step completely unchanged: -im*dz is already the
    correct generalization for a complex dz (real dt, used by
    evolve_and_measure_tdvp/quench_tdvp above, is just the dz.imag==0
    special case) -- Julia's `-im*dz` doesn't care whether dz is real or
    complex.

    Honors self.tevol_method exactly as the two real-time entry points
    above do, via the same _check_tevol_method: "TDVP_GSE" runs one-site
    TDVP with (when do_gse) a preceding global subspace expansion,
    mirroring tdz.py's own generic TDVP_GSE branch, and "TEBD"/"MPO"
    raise. Before this, tdz.py's julia_live branch returned here
    unconditionally and silently ran two-site TDVP whatever
    self.tevol_method said -- so on one and the same chain
    submode="TDZ" and evolve_and_measure() disagreed about which
    integrator "TDVP_GSE" meant, and "TEBD" raised for one and quietly
    ran TDVP for the other. self.tevol_method="AUTO" is accepted here too
    and reduces to plain "TDVP" (the else branch below already treats
    anything other than "TDVP_GSE" that way) -- TDZ has no TEBD
    counterpart to try first, see tdz.py's own docstring for why."""
    _check_tevol_method(self,methods=("TDVP","TDVP_GSE","AUTO"))
    from .juliasession import Main as Mainjl
    from .mps import MPS
    if self.tevol_method=="TDVP_GSE":
        psi = wf.jlmps
        if do_gse:
            psi = Mainjl.gse_expand(Hop.jlmpo,psi,
                    int(self.tdvp_gse_krylov_order),self.tdvp_gse_cutoff,
                    self.maxm)
        psi = Mainjl.tdvp_step_onesite(Hop.jlmpo,psi,dz,self.cutoff,self.maxm)
    else:
        psi = Mainjl.tdvp_step(Hop.jlmpo,wf.jlmps,dz,self.cutoff,self.maxm)
    return MPS(psi,MBO=self)
