# Non-Hermitian DMRG on the live Julia backend (itensor_version=
# "julia_live"), cross-checked against exact diagonalization on a
# Hatano-Nelson-style chain: *asymmetric* hopping (t_right != t_left)
# plus a staggered imaginary on-site potential, so the Hamiltonian is
# genuinely non-symmetric (H^T != H) as well as non-Hermitian.
#
# Why that matters, and why this example exists as its own script rather
# than another row in ../nhdmrg_VS_ED_VS_arnoldi: julia_live doesn't run
# one of this codebase's own ports of ITensorNHDMRG.jl, it calls the real
# package -- whose "adjoint" sweep is against swapprime(H,0=>1) == H^T,
# not H^dagger. So its left vector satisfies the *transpose* eigenvalue
# equation, while dmrgpy's convention throughout (nhdmrg.py) is the
# adjoint one, H^dagger|psil> = conj(lambda)|psil> with <psil|psir> = 1.
# mpsjulialive/nhdmrg.jl's nh_biorthogonal_pair reconciles the two with a
# complex conjugation. On a complex-*symmetric* Hamiltonian -- which
# every other non-Hermitian example in this folder happens to be, since
# their hoppings are symmetric and every non-Hermitian piece is diagonal
# -- the two conventions coincide up to exactly that conjugation, so the
# fix is invisible there and only a model like this one can show it.
#
# This model also exercises the second julia_live-specific fix in that
# file: the spectrum here is a complex-conjugate pair tying for the
# smallest real part, and nothing inside ITensorNHDMRG ties its left
# solve to whichever member the right solve picked, so the two used to
# converge to different eigenvalues (deterministically, on every attempt
# from every random start -- their overlap came out exactly 0).
# nhdmrg_solve() breaks the tie by re-solving against exp(i*theta)*H for
# a small theta, which leaves every eigenvector untouched and rotates the
# spectrum just enough to separate the pair's real parts.
#
# Sweeps the hopping asymmetry t_left (t_right fixed at 1), and plots the
# targeted eigenvalue (smallest real part, dmrgpy's non-Hermitian
# "ground state" convention) against ED, plus both eigen-residuals --
# the only meaningful convergence certificate here, since a non-Hermitian
# "energy" is not a variational bound.

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt

from dmrgpy import fermionchain

n = 4              # small enough that ED gives the full reference spectrum
t_right = 1.0
gamma = 0.7        # staggered imaginary potential
t_lefts = [1.0, 0.8, 0.6, 0.4, 0.2]


def build_chain(t_left, backend):
    fc = fermionchain.Fermionic_Chain(n)
    if backend == "julia": fc.setup_julia()
    else: fc.setup_python()
    h = 0
    for i in range(n-1): # asymmetric hopping: this is what makes H^T != H
        h = h + t_right*fc.Cdag[i]*fc.C[i+1] + t_left*fc.Cdag[i+1]*fc.C[i]
    for i in range(n): # staggered imaginary potential
        h = h + 1j*gamma*(-1)**i*fc.Cdag[i]*fc.C[i]
    for i in range(n-1): # interaction
        h = h + 0.5*(fc.N[i]-0.5)*(fc.N[i+1]-0.5)
    fc.set_hamiltonian(h)
    return fc,h


def norm(wf): return abs(wf.dot(wf))**0.5


lam_ed,lam_jl,res_r,res_l = [],[],[],[]

for t_left in t_lefts:
    # exact reference: the full non-Hermitian matrix, sorted by real part
    fc,h = build_chain(t_left,"python")
    es_ed = fc.get_excited(mode="ED",n=4)

    fc,h = build_chain(t_left,"julia")
    e,psil,psir = fc.nhdmrg() # (eigenvalue, left, right), <psil|psir>=1

    # right eigenvector: H|psir> = lambda|psir>
    rr = norm(h*psir - e*psir)
    # left eigenvector, in dmrgpy's ADJOINT convention (this is the check
    # that fails outright without nh_biorthogonal_pair's conjugation --
    # it sits at ~2 instead of ~1e-15)
    rl = norm(h.get_dagger()*psil - np.conj(e)*psil)

    # the returned eigenvalue must be one of ED's, at the smallest real
    # part (either member of a conjugate pair sharing it is acceptable)
    assert abs(e.real-es_ed[0].real) < 1e-6
    assert min(abs(e-x) for x in es_ed) < 1e-6
    assert abs(psil.dot(psir)-1.0) < 1e-8   # biorthonormal
    assert rr < 1e-3 and rl < 1e-3          # a genuine eigenpair, both sides

    lam_ed.append(es_ed[0]) ; lam_jl.append(e)
    res_r.append(rr) ; res_l.append(rl)
    print("t_left = %4.2f   lambda_ED = %12.8f%+12.8fj   lambda_julia = "
          "%12.8f%+12.8fj   res_R = %.2e   res_L = %.2e"
          %(t_left,es_ed[0].real,es_ed[0].imag,e.real,e.imag,rr,rl))

lam_ed = np.array(lam_ed) ; lam_jl = np.array(lam_jl)

fig,(ax_e,ax_r) = plt.subplots(1,2,figsize=(11,4.2))

ax_e.plot(t_lefts,lam_ed.real,"o-",label=r"Re $\lambda$, ED")
ax_e.plot(t_lefts,lam_jl.real,"x--",label=r"Re $\lambda$, NH-DMRG (julia)")
# ED and NH-DMRG can land on opposite members of a conjugate pair, so
# compare |Im| rather than Im -- both are equally valid answers here
ax_e.plot(t_lefts,np.abs(lam_ed.imag),"s-",label=r"|Im $\lambda$|, ED")
ax_e.plot(t_lefts,np.abs(lam_jl.imag),"+--",label=r"|Im $\lambda$|, NH-DMRG (julia)")
ax_e.set_xlabel(r"$t_\mathrm{left}$   ($t_\mathrm{right}=1$)")
ax_e.set_ylabel(r"$\lambda$")
ax_e.set_title("Smallest-real-part eigenvalue")
ax_e.legend(fontsize=8) ; ax_e.grid(alpha=0.3)

ax_r.semilogy(t_lefts,res_r,"o-",
        label=r"$\||H|\psi_R\rangle-\lambda|\psi_R\rangle\|$")
ax_r.semilogy(t_lefts,res_l,"s-",
        label=r"$\||H^\dagger|\psi_L\rangle-\bar\lambda|\psi_L\rangle\|$")
ax_r.set_xlabel(r"$t_\mathrm{left}$   ($t_\mathrm{right}=1$)")
ax_r.set_ylabel("eigen-residual")
ax_r.set_title("Biorthogonal pair: both residuals")
ax_r.legend(fontsize=8) ; ax_r.grid(alpha=0.3)

fig.tight_layout()
plt.show()
