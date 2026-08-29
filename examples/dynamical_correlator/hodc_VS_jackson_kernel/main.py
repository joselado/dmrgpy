# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# HODC vs Jackson: two reconstructions of the SAME KPM-DMRG moments.
#
# kernel="hodc" is the high-order delta-Chebyshev reconstruction of
# arXiv:2512.03149 (see algebra/kpm.py). Every classic KPM kernel --
# Jackson included -- damps the moments, which is the same as convolving
# with a positive kernel, and no positive kernel is better than
# second-order accurate: the error at a smooth point decays as O(N^-2) in
# the number of moments, full stop. HODC gives up positivity in exchange
# for an O(eta^m) regularized delta, and costs nothing extra to evaluate,
# because it consumes the same undamped moments.
#
# The claim is asymptotic and specifically about SMOOTH points, so this
# script is built around a case where "smooth" is checkable rather than
# asserted:
#
#   * The model is the XX chain, H = sum_i (Sx_i Sx_j + Sy_i Sy_j), which
#     Jordan-Wigner maps to free fermions. <Sz_i delta(w-H+E0) Sz_i> is
#     then a particle-hole spectrum whose lines and weights follow from
#     one LxL single-particle diagonalization -- so the exact answer for
#     this same finite chain, and the L->infinity continuum it converges
#     to, are both a few lines of numpy. No ED, no finite-size mismatch
#     that has to be argued away.
#
#   * That continuum is smooth except at a van Hove singularity at w=1
#     and at the band edges w=0,2. The error is therefore measured twice:
#     once on a window that avoids the singularity (where the theory
#     predicts a gain) and once on a window that contains it (where it
#     predicts much less). Reporting only the first would be cheating.
#
# The DMRG side is deliberately modest -- kpmmaxm=32 -- because that is
# the regime the question actually arises in, and it is where the answer
# stops being the paper's answer. Each kernel is therefore run on two
# moment sources: the DMRG moments, and the exact ones (free fermions
# again -- a matter of milliseconds). On exact moments HODC behaves as
# advertised. On the DMRG moments it eventually loses to Jackson, and the
# reason is visible in panel 3: MPS truncation error in the Chebyshev
# recursion grows exponentially in k, and Jackson's damping factors --
# which go to zero exactly at k=N -- happen to suppress precisely the
# moments that are worst. HODC consumes the moments undamped, so it
# inherits that noise; its own e^{-k*eta} tail decay is the only thing
# holding it back, which is why raising eta buys noise immunity (panel 4)
# and why the useful N is capped by moment accuracy, not by patience.
#
# Measured alongside this script, on the same chain: kpm_accelerate does
# not cause that horizon. With it on the moments pass 1e-3 error at k=38
# and end up O(1e2) by k=1200; with it off they pass 1e-3 at k=39 and
# stay below 1. So the accelerated recursion only makes the late blow-up
# violent -- the onset is set by kpmmaxm, and turning kpm_accelerate off
# buys stability, not accuracy.
#
# One moment run feeds everything: prefixes of a single length-PMAX
# moment array are exactly the moments a shorter run would have produced,
# so the error-vs-N curves compare kernels and not two DMRG runs that
# converged slightly differently.
import time
import numpy as np
import matplotlib.pyplot as plt

from dmrgpy import spinchain
from dmrgpy import kpmdmrg

L = 48          # chain length (even -> no zero mode -> unique ground state)
KPMMAXM = 32    # bond dimension of the Chebyshev vectors: modest on purpose
PMAX = 1200     # moments computed once; every N below is a prefix of these
SITE = L//2-1   # middle site
KPM_SCALE = 0.9 # margin left inside [-1,1] by the spectrum rescaling.
# Raised from dmrgpy's 0.7 default deliberately, and not because the band
# edges are poor -- they are exact to ~1e-5 here, since the lower one is
# just the (well converged) ground state and the XX spectrum is symmetric.
# It is the Chebyshev recursion itself that is only conditionally stable
# under MPS truncation: at kpmmaxm=32 and a thousand-plus moments, the
# truncation error each step introduces gets amplified by T_k outside
# [-1,1], and at kpm_scale=0.7 (spectrum filling |x|<=0.71) that
# eventually trips the backend's own divergence guard. More margin damps
# the amplification. The cost is resolution: the reconstruction's width
# at fixed N is proportional to kpm_scale -- which applies identically to
# both kernels, so the comparison below is unaffected.


# ----------------------------------------------------------------------
# exact free-fermion reference
# ----------------------------------------------------------------------
def sticks(n):
    """Exact <Sz_i delta(w-H+E0) Sz_i> of the open XX chain of n sites at
    half filling, as a list of (energy, weight) lines.

    Sz_i = c^dag_i c_i - 1/2 with no Jordan-Wigner string, so acting on
    the Slater-determinant ground state it produces exactly the
    particle-hole pairs |p,h>, with energy eps_p-eps_h and weight
    |phi_p(i)|^2 |phi_h(i)|^2."""
    hop = np.zeros((n,n))
    for i in range(n-1): hop[i,i+1] = hop[i+1,i] = 0.5 # Sx Sx + Sy Sy
    e,v = np.linalg.eigh(hop)
    occ = e<0.0                       # fill the negative-energy modes
    phi2 = v[n//2-1,:]**2
    om = (e[~occ][None,:]-e[occ][:,None]).ravel()
    wt = (phi2[occ][:,None]*phi2[~occ][None,:]).ravel()
    return om,wt,np.sum(e[occ])


def smooth_density(om,wt,ws,sigma):
    """Lines -> a density on the grid ws, Gaussian-broadened by sigma.
    Done by histogram + convolution so a million lines stay cheap."""
    dw = ws[1]-ws[0]
    edges = np.concatenate([ws-dw/2.,[ws[-1]+dw/2.]])
    h,_ = np.histogram(om,bins=edges,weights=wt)
    t = np.arange(-int(6*sigma/dw),int(6*sigma/dw)+1)*dw
    g = np.exp(-t**2/(2*sigma**2)) ; g /= np.sum(g)*dw # int g dw = 1
    return np.convolve(h,g,mode="same") # sum_j wt_j G(w-om_j)


# ----------------------------------------------------------------------
# the DMRG side: one moment run, reused by every kernel below
# ----------------------------------------------------------------------
sc = spinchain.Spin_Chain(["S=1/2"]*L)
sc.itensor_version = "python" ; sc.setup_python()
h = 0
for i in range(L-1): h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
sc.set_hamiltonian(h)
sc.maxm,sc.nsweeps,sc.kpmmaxm = 60,8,KPMMAXM
sc.kpm_scale = KPM_SCALE

om_L,wt_L,e0_exact = sticks(L)
t0 = time.time()
e0 = sc.gs_energy()
print("ground state   DMRG %.10f   exact %.10f   (%.0f s)"%(
    e0,e0_exact,time.time()-t0))

# delta is what sets the moment count: N = (emax-emin)/delta * kpm_n_scale
delta = 2*abs(e0_exact)*sc.kpm_n_scale/PMAX
t0 = time.time()
mus,emin,emax,scale,npol,delta = sc.get_dynamical_correlator_moments(
        name=(sc.Sz[SITE],sc.Sz[SITE]),delta=delta)
print("%d KPM moments at kpmmaxm=%d in %.0f s (%.2f s/moment)"%(
    npol,KPMMAXM,time.time()-t0,(time.time()-t0)/npol))

# the same moments, computed exactly from the free-fermion lines
xj = (emin+om_L-(emin+emax)/2.)*scale
mus_exact = np.array([np.sum(wt_L*np.cos(k*np.arccos(xj)))
                      for k in range(npol)])
print("moment sum rule: mu_0 = %.8f (exact <Sz^2> = 0.25)"%mus[0].real)


# ----------------------------------------------------------------------
# reconstruct prefixes of those moments with both kernels
#
# Two moment SOURCES are compared, not one: the DMRG moments just
# computed, and the exact ones. That is what separates the kernel
# question ("is a high-order delta better than a damped one?") from the
# MPS question ("are these moments good enough to tell?"), and the second
# turns out to matter at least as much at kpmmaxm=32.
# ----------------------------------------------------------------------
WINDOWS = {"smooth (no van Hove)":(0.30,0.85),
           "with van Hove at w=1":(0.30,1.70)}
SOURCES = {"exact":mus_exact,"DMRG":mus.real}
NS = [n for n in (150,300,600,1200) if n<=npol]

# reference: the L->infinity continuum this finite chain converges to.
# At the resolutions probed here the L=48 lines are already within
# ~2e-3 of it, well below every error compared below.
ws_ref = np.linspace(0.0,2.05,2100)
om_inf,wt_inf,_ = sticks(1024)
# (the line weights already sum to <Sz^2>=1/4 for every L, so no
# renormalization is needed when swapping L=48 for L=1024 here)
rho_inf = smooth_density(om_inf,wt_inf,ws_ref,0.004)


def reconstruct(src,kern,n,eta=None,grid=None):
    _,y = kpmdmrg.dynamical_correlator_from_moments(
            SOURCES[src][:n],emin,emax,scale,n,
            ws_ref if grid is None else grid,kernel=kern,
            delta=(emax-emin)*sc.kpm_n_scale/n,hodc_eta=eta)
    return y.real

def rel_error(y,window):
    lo,hi = WINDOWS[window]
    s = (ws_ref>lo)&(ws_ref<hi)
    return np.linalg.norm((y-rho_inf)[s])/np.linalg.norm(rho_inf[s])

errors,spectra = {},{}
for src in SOURCES:
    for kern in ("jackson","hodc"):
        for n in NS:
            y = reconstruct(src,kern,n)
            spectra[(src,kern,n)] = y
            for w in WINDOWS:
                errors.setdefault((w,src,kern),[]).append(rel_error(y,w))

for w in WINDOWS:
    print("\n%s   relative L2 error vs the exact continuum"%w)
    print("  %6s %11s %11s %7s   %11s %11s %7s"%(
        "N","exact:jack","exact:hodc","gain","DMRG:jack","DMRG:hodc","gain"))
    for i,n in enumerate(NS):
        e = {k:errors[(w,)+k][i] for k in
             [("exact","jackson"),("exact","hodc"),
              ("DMRG","jackson"),("DMRG","hodc")]}
        print("  %6d %11.3e %11.3e %6.1fx   %11.3e %11.3e %6.1fx"%(
            n,e[("exact","jackson")],e[("exact","hodc")],
            e[("exact","jackson")]/e[("exact","hodc")],
            e[("DMRG","jackson")],e[("DMRG","hodc")],
            e[("DMRG","jackson")]/e[("DMRG","hodc")]))

nbig = NS[-1]
print()
# Normalization: int S^zz(w) dw must be <Sz^2> = 1/4. Integrated on a
# grid wider than the plotted one, because at the coarsest N here the
# reconstruction is ~0.5 wide and a [0,2.05] window would lose several
# percent of the weight over the band edges -- i.e. it would be measuring
# the window rather than the spectrum.
#
# On the exact moments this is a pure plumbing check (scale, shift, units)
# and is asserted. On the DMRG moments it is a diagnostic instead, and a
# revealing one: it is *not* protected from the moment noise. Each kernel
# integrates to one, so the total weight is fixed by mu_0 -- but only in
# exact arithmetic; what actually gets integrated is sum_k mu_k * int
# nu_k(E) dE, and a residual of ~1e-4 in the k-th term is enough to move
# the answer by 10% once |mu_k| has grown to O(100). Watch the HODC
# column walk away from 1/4 as N passes the moment horizon.
ws_wide = np.linspace(-2.0,4.0,4000)
for kern in ("jackson","hodc"):
    w = np.trapezoid(reconstruct("exact",kern,nbig,grid=ws_wide),ws_wide)
    print("total weight, exact moments, N=%d, %-8s %.6f (exact 0.25)"%(
        nbig,kern,w))
    assert abs(w-0.25) < 5e-3
for kern in ("jackson","hodc"):
    ws_ = [np.trapezoid(reconstruct("DMRG",kern,n,grid=ws_wide),ws_wide)
           for n in NS]
    print("total weight, DMRG  moments, %-8s "%kern
          +"  ".join("N=%d %.5f"%(n,w) for n,w in zip(NS,ws_)))

# The kernel claim itself, on clean moments and a smooth window: this is
# a property of HODC rather than of this DMRG run, so it is the part
# worth asserting. It is checked as a best-over-N statement because both
# kernels eventually hit the same L=48 finite-size floor (~5e-3 against
# an L->infinity reference), which caps the asymptotic gain from above
# and would make a fixed-N threshold brittle. The unclipped, purely
# numerical version of this claim -- HODC an order of magnitude better on
# a genuinely smooth density -- lives in tests/test_hodc_kernel.py.
w0 = "smooth (no van Hove)"
gains = [ej/eh for ej,eh in zip(errors[(w0,"exact","jackson")],
                               errors[(w0,"exact","hodc")])]
print("exact-moment gain, smooth window: "
      +", ".join("N=%d %.1fx"%(n,g) for n,g in zip(NS,gains)))
assert max(gains) > 2.0, "HODC should beat Jackson on a smooth continuum"


# ----------------------------------------------------------------------
# HODC's one free knob, eta, swept at fixed N -- the practical question
# once the moments are known to be imperfect. nu_k ~ q^k decays like
# exp(-k*eta), so eta is simultaneously the resolution and the only
# damping HODC applies to the high-k moments.
# ----------------------------------------------------------------------
# (the dmrgpy default is spliced into the grid so the sweep's own
# minimum and the default-eta number below are directly comparable)
netas = np.sort(np.append(np.linspace(1.0,24.0,35),
                          sc.kpm_n_scale/sc.kpm_scale))
eta_scan = {src:[rel_error(reconstruct(src,"hodc",nbig,eta=ne/nbig/scale),w0)
                 for ne in netas] for src in SOURCES}
print()
for src in SOURCES:
    i = int(np.argmin(eta_scan[src]))
    ej = errors[(w0,src,"jackson")][-1]
    print("%-5s moments, N=%d: HODC is best at N*eta=%.1f (err %.3e, "
          "%.1fx Jackson); at the default N*eta=%.1f it is %.3e"%(
          src,nbig,netas[i],eta_scan[src][i],ej/eta_scan[src][i],
          sc.kpm_n_scale/sc.kpm_scale,
          errors[(w0,src,"hodc")][-1]))


# ----------------------------------------------------------------------
# plots
# ----------------------------------------------------------------------
fig,ax = plt.subplots(2,2,figsize=(12.5,9))

a = ax[0,0]
nshow = NS[1]
a.plot(ws_ref,rho_inf,color="k",lw=2.4,label="exact ($L\\to\\infty$)")
a.plot(ws_ref,spectra[("DMRG","jackson",nshow)],color="C0",lw=1.2,
       label="Jackson, N=%d"%nshow)
a.plot(ws_ref,spectra[("DMRG","hodc",nshow)],color="C3",lw=1.2,
       label="HODC m=6, N=%d"%nshow)
a.axvline(1.0,color="gray",ls=":",lw=1)
a.text(1.03,0.015,"van Hove",color="gray",fontsize=8)
lo,hi = WINDOWS["smooth (no van Hove)"]
a.axvspan(lo,hi,color="C2",alpha=0.10)
a.set_xlabel("$\\omega$") ; a.set_ylabel("$S^{zz}_{ii}(\\omega)$")
a.set_xlim(0,2.05) ; a.legend(fontsize=8)
a.set_title("XX chain, L=%d, kpmmaxm=%d, DMRG moments\n"
            "(shaded: the smooth error window)"%(L,KPMMAXM),fontsize=9)

a = ax[0,1]
for src,ls,mk in (("exact","-","o"),("DMRG","--","s")):
    for kern,col in (("jackson","C0"),("hodc","C3")):
        a.loglog(NS,errors[(w0,src,kern)],ls+mk,color=col,ms=4,
                 label="%s moments, %s"%(src,kern))
nn = np.array(NS,dtype=float)
a.loglog(nn,errors[(w0,"exact","jackson")][0]*(nn/nn[0])**-2.,":",
         color="gray",lw=1,label="$N^{-2}$")
a.set_xlabel("number of moments $N$")
a.set_ylabel("relative $L_2$ error, smooth window")
a.legend(fontsize=7)
a.set_title("solid = kernel alone; dashed = same kernel on\n"
            "moments an MPS at kpmmaxm=%d actually produces"%KPMMAXM,
            fontsize=9)

a = ax[1,0]
kk = np.arange(npol)
a.semilogy(kk,np.abs(mus[:npol].real-mus_exact)+1e-18,color="C4",lw=0.7,
           label="$|\\mu_k^{\\rm DMRG}-\\mu_k^{\\rm exact}|$")
a.axhline(abs(mus_exact[0]),color="k",ls=":",lw=1,
          label="$\\mu_0$ (the moments' own scale)")
for n,c in zip(NS,("C7","C8","C9","C1")): a.axvline(n,color=c,lw=0.8,ls="--")
a.set_xlabel("moment index $k$")
a.set_ylabel("absolute moment error")
a.legend(fontsize=7)
a.set_title("why the dashed curves turn up: MPS truncation noise\n"
            "grows exponentially in $k$ (dashed lines: the $N$ above)",
            fontsize=9)

a = ax[1,1]
for src,ls,mk in (("exact","-","o"),("DMRG","--","s")):
    a.semilogy(netas,eta_scan[src],ls+mk,color="C3",ms=3,
               label="HODC m=6, %s moments"%src)
    a.axhline(errors[(w0,src,"jackson")][-1],color="C0",ls=ls,lw=1,
              label="Jackson, %s moments"%src)
a.axvline(sc.kpm_n_scale/sc.kpm_scale,color="gray",ls=":",
          label="dmrgpy default $N\\eta$")
a.set_xlabel("$N\\eta$   (regularization width $\\times$ moment count)")
a.set_ylabel("relative $L_2$ error, smooth window")
a.legend(fontsize=7)
a.set_title("HODC's one knob at N=%d. $\\nu_k\\propto q^k$ decays like\n"
            "$e^{-k\\eta}$, so a larger $\\eta$ also buys noise damping"%nbig,
            fontsize=9)

fig.tight_layout()
fig.savefig("hodc_VS_jackson_kernel.png",dpi=130)
plt.show()
