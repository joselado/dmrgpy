"""Figure 6, excited-states_heff: the exchange J of the Hubbard dimer,
fitted by get_heff in the four lowest states, against the exact dimer
splitting and against 4t^2/U.

Run with  bash docs/figures/run.sh excited-states_heff

The dimer is H = t sum_s (c^dag_0s c_1s + h.c.) + U sum_i (n_iup-1/2)(n_idn-1/2)
+ lam (N-2)^2. The last term pins half filling: without it the N=1 and
N=3 states fall below the triplet for U < 2t, and the four lowest states
are not the singlet and the triplet. It commutes with H and with the spin
operators, so it does not change the fit inside the N=2 manifold.

What the fit returns. get_heff fits P H P = c + J P S0.S1 P with the
projected operators. The triplet is purely covalent, so <S0.S1> = 1/4
there, but the singlet is a mix of the covalent singlet (weight w) and
the symmetric doubly occupied state, so <S0.S1> = -3w/4 there. The fit is
then J = (E_T - E_S)/(1/4 + 3w/4), which is larger than the splitting
E_T - E_S and joins it as w -> 1 at large U. With the 2x2 problem
[[0, -2t], [-2t, U]] in the singlet sector, the amplitude ratio of the
ionic to the covalent part is b/a = -E_S/(2t), with
E_S = (U - sqrt(U^2+16t^2))/2, and w = 1/(1+(b/a)^2). This analytic
j_proj is the check of the code; j_exact and j_strong are the physics.

Arrays in data.npz:
    U         (k,)    Hubbard U, in units of t
    jfit      (k,3)   fitted couplings of SxSx, SySy, SzSz
    j_exact   (k,)    exact dimer splitting (sqrt(U^2+16t^2)-U)/2
    j_strong  (k,)    strong-coupling exchange 4t^2/U
    j_proj    (k,)    exact value of the projected fit, j_exact/(1/4+3w/4)
    w         (k,)    weight of the singly occupied (covalent) part of the
                      singlet ground state
    params            JSON string with every parameter (see common.py)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common
from dmrgpy import fermionchain

params = common.start(dict(
    model="Hubbard dimer, Spinful_Fermionic_Chain(2), PH-symmetric U term,"
          " half filling pinned by lam*(N-2)^2",
    t=1.0,
    U=[float(u) for u in np.geomspace(1, 20, 25)],
    lam=20.0,
    n_states=4,
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
jfit = []
for U in us:
    fc, ops = dimer(U)
    j = np.real(np.array(fc.get_heff(operators=ops, n=params["n_states"],
                                     mode=params["mode"])))
    jfit.append(j)
    je = 0.5*(np.sqrt(U**2 + 16*t**2) - U)
    print("U=%6.3f  J_xx=%.10f J_yy=%.10f J_zz=%.10f  exact=%.10f  "
          "4t^2/U=%.6f"
          % (U, j[0], j[1], j[2], je, 4*t**2/U),
          flush=True)
jfit = np.array(jfit)
j_exact = 0.5*(np.sqrt(us**2 + 16*t**2) - us)
j_strong = 4*t**2/us
e_s = 0.5*(us - np.sqrt(us**2 + 16*t**2))
w = 1.0/(1.0 + (e_s/(2*t))**2)
j_proj = j_exact/(0.25 + 0.75*w)
print("max over U of |J_fit - J_proj| = %.3e  (the check)"
      % np.max(np.abs(jfit - j_proj[:, None])))
print("spread of the three couplings, max over U = %.3e"
      % np.max(jfit.max(axis=1) - jfit.min(axis=1)))
print("relative J_fit/J_exact - 1 at U=%g: %.3e, at U=%g: %.3e"
      % (us[0], jfit[0, 2]/j_exact[0]-1, us[-1], jfit[-1, 2]/j_exact[-1]-1))

common.save(__file__, params, U=us, jfit=jfit, j_exact=j_exact,
            j_strong=j_strong, j_proj=j_proj, w=w)
