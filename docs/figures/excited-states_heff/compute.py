"""Figure 6, excited-states_heff: the exchange J of the Hubbard dimer,
fitted by get_heff in the four lowest states, against the exact dimer
splitting E_T - E_S and against its strong-coupling limit 4t^2/U.

Run with  bash docs/figures/run.sh excited-states_heff

The dimer is H = t sum_s (c^dag_0s c_1s + h.c.) + U sum_i (n_iup-1/2)(n_idn-1/2)
+ lam (N-2)^2. The last term pins half filling: without it the N=1 and
N=3 states fall below the triplet for U < 2t, and the four lowest states
are not the singlet and the triplet. It commutes with H and with the spin
operators, so it does not change the fit inside the N=2 manifold.

Two fits. With projector=single_occupancy_projector(fc), the four
eigenstates are projected on the singly occupied states and Lowdin
orthonormalized, and get_heff fits H_eff = sum_n E_n |n~><n~|, so the
fitted J is exactly E_T - E_S (jfit). Without the projector, the spin
operators are fitted inside the four eigenstates themselves; the singlet
is a mix of the covalent singlet (weight w) and the symmetric doubly
occupied state, where S0.S1 = 0, so that fit gives
J = (E_T - E_S)/(1/4 + 3w/4), above the splitting at small U (jfit_noproj).
With the 2x2 problem [[0, -2t], [-2t, U]] in the singlet sector,
E_S = (U - sqrt(U^2+16t^2))/2 and w = 1/(1+(E_S/2t)^2).

Arrays in data.npz:
    U            (k,)    Hubbard U, in units of t
    jfit         (k,3)   fitted couplings of SxSx, SySy, SzSz, with the
                         singly occupied projector
    jfit_noproj  (k,3)   the same without the projector
    j_exact      (k,)    exact dimer splitting (sqrt(U^2+16t^2)-U)/2
    j_strong     (k,)    strong-coupling exchange 4t^2/U
    j_proj       (k,)    exact value of the fit without projector,
                         j_exact/(1/4+3w/4)
    w            (k,)    weight of the singly occupied (covalent) part of
                         the singlet ground state
    params               JSON string with every parameter (see common.py)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import effectivehamiltonian, fermionchain

params = common.start(dict(
    model="Hubbard dimer, Spinful_Fermionic_Chain(2), PH-symmetric U term,"
          " half filling pinned by lam*(N-2)^2",
    t=1.0,
    U=[float(u) for u in np.geomspace(1, 20, 25)],
    lam=20.0,
    n_states=4,
    projector="single_occupancy_projector",
    mode="ED",
))
t, lam = params["t"], params["lam"]


def dimer(U):
    fc = fermionchain.Spinful_Fermionic_Chain(2)
    h = 0
    h = h + t*fc.Cdagup[0]*fc.Cup[1] + t*fc.Cdagdn[0]*fc.Cdn[1]
    h = h + h.get_dagger()
    for i in range(2):
        h = h + U*(fc.Nup[i]-0.5)*(fc.Ndn[i]-0.5)
    ntot = fc.Ntot[0] + fc.Ntot[1]
    h = h + lam*(ntot-2)*(ntot-2)
    fc.set_hamiltonian(h)
    ops = [fc.Sx[0]*fc.Sx[1], fc.Sy[0]*fc.Sy[1], fc.Sz[0]*fc.Sz[1]]
    return fc, ops


us = np.array(params["U"])
jfit, jfit_noproj = [], []
for U in us:
    fc, ops = dimer(U)
    proj = effectivehamiltonian.single_occupancy_projector(fc)
    j = np.real(np.array(fc.get_heff(operators=ops, n=params["n_states"],
                                     mode=params["mode"], projector=proj)))
    j0 = np.real(np.array(fc.get_heff(operators=ops, n=params["n_states"],
                                      mode=params["mode"])))
    jfit.append(j)
    jfit_noproj.append(j0)
    je = 0.5*(np.sqrt(U**2 + 16*t**2) - U)
    print("U=%6.3f  J_xx=%.10f J_yy=%.10f J_zz=%.10f  no projector=%.10f  "
          "exact=%.10f  4t^2/U=%.6f"
          % (U, j[0], j[1], j[2], j0[2], je, 4*t**2/U),
          flush=True)
jfit = np.array(jfit)
jfit_noproj = np.array(jfit_noproj)
j_exact = 0.5*(np.sqrt(us**2 + 16*t**2) - us)
j_strong = 4*t**2/us
e_s = 0.5*(us - np.sqrt(us**2 + 16*t**2))
w = 1.0/(1.0 + (e_s/(2*t))**2)
j_proj = j_exact/(0.25 + 0.75*w)
print("max over U of |J_fit - (E_T-E_S)| = %.3e  (the check)"
      % np.max(np.abs(jfit - j_exact[:, None])))
print("max over U of |J_fit(no projector) - J_proj| = %.3e"
      % np.max(np.abs(jfit_noproj - j_proj[:, None])))
print("spread of the three couplings, max over U = %.3e"
      % np.max(jfit.max(axis=1) - jfit.min(axis=1)))
print("without projector, J/(E_T-E_S) - 1 at U=%g: %.3e, at U=%g: %.3e"
      % (us[0], jfit_noproj[0, 2]/j_exact[0]-1,
         us[-1], jfit_noproj[-1, 2]/j_exact[-1]-1))

common.save(__file__, params, U=us, jfit=jfit,
            jfit_noproj=jfit_noproj, j_exact=j_exact,
            j_strong=j_strong, j_proj=j_proj, w=w)
