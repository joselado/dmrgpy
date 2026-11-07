"""A ground state of definite fermion parity for the STM spectrum of a
magnetic impurity on a superconductor.

Every tunneling vertex X of ysr.tunneling_operators adds or removes one
electron, so it changes the total fermion parity, and on a ground state of
definite parity <GS|X|GS> = 0 exactly. The unrestricted DMRG solve that the
TASK correlator starts from does not know this. Near the singlet-doublet
transition the two parities are nearly degenerate, and a solve from a
random start can stop at a mixture a|even> + b|odd>. On that state the
elastic weight |<X^dag><X>| is nonzero, the subgap pole splits and moves,
and ysrcurrent.rate_parameters halves u2. stm_poles and rate_parameters
warn about it, but a warning does not give the right spectrum.

The cure is to solve inside one parity sector. The impurity sites each
hold exactly one electron (the penalty Upen of ysr.kondo_chain), so total
Sz fixes the parity: for S=1/2, 2Sz = 0 is even total parity (the
screened singlet, odd bath) and 2Sz = 1 odd (the free doublet). For a spin
S held as 2S sites the bath parity is (2Sz - 2S) mod 2, and with Dz or B
the lowest state of a parity can sit at any |2Sz| <= 2S+1. For N
surrogate spins 2Sz = N + N_bath mod 2. parity_ground_state solves the
lowest state of each 2Sz sector on a copy of the chain, with the warm
starts that keep the impurity charges right (ysr._seed_impurity_charge on
v3 for S > 1/2, surrogate._warm_sector_solve for several spins) and the
charge checks after the solve. It keeps the copy of the lowest sector and
calls promote_to_dense() on it: the state leaves the sector exactly, on
dense indices, and the chain hands it to every correlator
(groundstate.ground_state_on_session), so a charge-changing X can act on
it. ITensor Index identity is per session, so the state is used on the
copy that produced it; the copy is what is returned.

With mode="ED" nothing of this is needed. ED diagonalizes the full
Hamiltonian, so its ground state is an eigenstate of the parity unless the
two parities are exactly degenerate, a fine-tuned point. And an ED chain
with a sector set restricts every operator to that sector, which a
charge-changing X leaves. So on mode="ED" the sector energies are computed
the same way, and the returned chain is a plain copy of fc.

What was checked (tests/test_ysr_parity.py, runs on d22, BLAS pinned to
one thread). TASK works on the promoted state on both backends: on the
4-site chain (J=0.6, Delta=0.1, U = 0 and 0.3, maxm=60, nsweeps=15) the
spectrum agrees with ED (CVM) to 2.6e-11 (v3) and 8.2e-10 ("python"), with
an elastic weight of 0 (v3) or 2e-29 and a missing weight below 7e-16.
Near the transition (8 sites, Lambda=4, J = 0.6 A_Lambda, Delta=0.1,
U=0.3, maxm=60, nsweeps=15, v3) four plain solves from fresh chains all
stopped at a mixture, with elastic weights of 1.2e-2, 9.7e-5, 3.8e-3 and
3.6e-6 of the total, subgap poles split by up to 1e-3 Delta and placed
down to 0.0169 Delta, and u2 between 0.153 and 0.33378. In the parity
sector the elastic weight is 0, there is one pole per bias sign at
0.0190312 Delta (ED: 0.01903116 Delta), and rate_parameters gives
u2 = 0.333787 and v2 = 0.196261 with no warning, the values of a plain
solve that happened to converge at nsweeps=40. The two sector solves cost
about 9 s there, against 7 to 11 s for one plain solve. On a surrogate
dimer next to a parity crossing (L=3, J=0.671, Delta=0.1, first odd state
1.709e-3 above the ground state by ED), three plain solves at maxm=64,
nsweeps=20 were 9.4e-5 to 5.2e-4 above the ED energy on v3 and their dI/dV
off by 1.1e-2 to 4.0e-2 of the peak (7e-5 to 8.8e-4 on "python"); the
parity state is 3.4e-12 (v3) and 2.3e-11 ("python") of the peak from ED.
That dimer has J12 = t_b = 0, so its two spins do not interact. With
J12 = t_b = 0.01 (maxm=64, nsweeps=20) the parity route is off ED by
2.5e-6 of the peak on "python" and 4.5e-9 on v3, where one plain v3 run
was 2.9e-5 off; the plain "python" route is 5.8e-7 off. The "python"
sector solve used to be the slow one here, 1.06e-4 of the peak at
nsweeps=20: its local solves stopped relative to the scale of H plus the
charge penalty, about 2e3 times that of H, so the state kept an
admixture of about 1e-4 in amplitude (1.4e-8 in weight) of the first
excited state of its sector (see pyitensor/dmrg.py's dmrg(penalty=)).

One thing the sector state changes besides the parity: TASK works in the
one-site tangent space of the state it gets, and that space is set by the
bond dimensions. A ground state of definite 2Sz has a smaller Schmidt rank
than the mixture of a degenerate pair +-2Sz that a plain solve returns.
For S=1, Dz=-0.05 and S=3/2, Dz=-0.04 (n=4 and 3, ground state at
2Sz = +-2) the tangent space of the sector state at the chain's cutoff
1e-12 was 562 and 205, against 786 and 297 for the plain state, and its
spectrum differed from ED by 8.7e-4 and 4.8e-2 of the peak at high
energy, while the subgap pole was exact (4e-10). With fc.cutoff = 0 the
tangent space fills the Hilbert space and the spectrum is ED's again
(2.6e-7 and 8e-11 on v3, 1.6e-7 and 2e-10 on "python"). This is the
tangent-space approximation of TASK, which a short chain hides by
saturating it; on a long chain both states are limited by maxm. The
chain's cutoff is used as it is; set fc.cutoff = 0 before the call to
saturate a short chain.
"""

import copy
import fractions

import numpy as np


def _is_surrogate(fc):
    return hasattr(fc, "impurities")


def default_sectors(fc):
    """The 2Sz sectors parity_ground_state solves by default, the ones of
    ysr.subgap_energies for a kondo_chain and of surrogate.sector_energies
    for a surrogate chain."""
    B = getattr(fc, "B_imp", 0.0)
    if _is_surrogate(fc):
        N = len(fc.impurities)
        return [0, 1] if B == 0.0 else list(range(-(N + 1), N + 2))
    ni = int(2*getattr(fc, "spin", fractions.Fraction(1, 2)))
    if getattr(fc, "Dz", 0.0) == 0.0 and B == 0.0:
        return [0, 1]
    sectors = list(range(0, ni + 2))
    if B != 0.0:
        sectors = [-s for s in sectors[::-1] if s > 0] + sectors
    return sectors


def _solve_sector(fc, sz, mode):
    """(copy of fc solved in the sector 2Sz = sz, its energy, whether it
    was solved by DMRG)."""
    c = copy.deepcopy(fc)
    c.set_conserved_sector(Sz=sz)
    dmrg = c.get_mode(mode=mode) == "DMRG"
    if _is_surrogate(fc):
        from . import surrogate
        if dmrg:
            e0 = surrogate._warm_sector_solve(c)
            q = surrogate._impurity_charge_error(c, mode)
            if q > 1e-6:
                raise RuntimeError("parity_ground_state: sector 2Sz=%d, the "
                    "DMRG state has sum_i <(N_i-1)^2> = %.2e on the spin "
                    "sites; raise nsweeps or maxm" % (sz, q))
        else:
            e0 = c.gs_energy(mode=mode)
        return c, float(np.real(e0)), dmrg
    from .. import ysr
    ni = int(2*getattr(fc, "spin", fractions.Fraction(1, 2)))
    if ni > 1 and dmrg and c.itensor_version == 3:
        ysr._seed_impurity_charge(c, ni)
    e0 = c.gs_energy(mode=mode)
    if ni > 1:
        ysr._check_impurity_charge(c, ni, mode)
    return c, float(np.real(e0)), dmrg


def parity_ground_state(fc, sectors=None, mode="DMRG"):
    """Ground state of fc solved inside its fermion-parity sector, on a
    copy of fc that is returned with the state promoted to dense indices,
    so that a correlator of a charge-changing operator (the tunneling
    vertex) starts from it. fc itself is left untouched.

    fc is a ysr.kondo_chain or a surrogate.surrogate_kondo_chain.
    sectors are the 2Sz values solved (default_sectors: [0, 1] at Dz = B =
    0, otherwise every 2Sz that can hold the lowest state of a parity).
    Each sector is solved on its own copy, the lowest is kept (on a
    degenerate multiplet the smallest |2Sz| within 1e-10 max(1, |E0|), the
    rule of ysr.subgap_energies; any member of a degenerate pair gives the
    same spin-summed spectrum), and its copy gets promote_to_dense().

    Returns (chain, info). info holds "sectors" (2Sz -> energy), "E0",
    "Sz_ground", "E_even" and "E_odd" (lowest energies of each total
    parity), "parity" ("even" or "odd", of the returned state), and for a
    kondo_chain "E_free", "E_screened" and "E_ysr" = E_free - E_screened as
    ysr.subgap_energies defines them; for a surrogate chain "E_bath_even"
    and "E_bath_odd". With mode="ED" (or a chain that the dispatch sends
    to ED) the energies are the ED sector energies and the chain is a
    plain copy of fc, whose ED ground state already has a definite parity
    (module docstring). Raises at Delta = 0, where an Sz sector does not
    fix the particle number, and on a backend with no sectors
    (itensor_version 2, "julia_live")."""
    if getattr(fc, "Delta_sc", 0.0) == 0.0:
        raise ValueError("parity_ground_state needs a gap, Delta > 0: at "
            "Delta = 0 the particle number is conserved too, so an Sz "
            "sector alone does not fix the state, and a solve stays at the "
            "particle number it started from")
    if sectors is None:
        sectors = default_sectors(fc)
    sectors = sorted(set(int(s) for s in sectors))
    if not any(s % 2 == 0 for s in sectors) or \
       not any(s % 2 == 1 for s in sectors):
        raise ValueError("sectors must contain both parities of 2Sz")
    ens, chains, dmrg = {}, {}, False
    for sz in sectors:
        c, e0, dmrg = _solve_sector(fc, sz, mode)
        ens[sz], chains[sz] = e0, c
    e0 = min(ens.values())
    tie = 1e-10*max(1.0, abs(e0)) # a degenerate multiplet: smallest |2Sz|
    szg = min((s for s, e in ens.items() if e <= e0 + tie),
              key=lambda s: (abs(s), s))
    info = dict(sectors=ens, E0=ens[szg], Sz_ground=szg,
                parity="even" if szg % 2 == 0 else "odd")
    def lowest(cond):
        v = [e for s, e in ens.items() if cond(s)]
        return min(v) if v else None
    info["E_even"] = lowest(lambda s: s % 2 == 0)
    info["E_odd"] = lowest(lambda s: s % 2 == 1)
    if _is_surrogate(fc):
        N = len(fc.impurities)
        info["E_bath_even"] = lowest(lambda s: (s - N) % 2 == 0)
        info["E_bath_odd"] = lowest(lambda s: (s - N) % 2 == 1)
    else:
        ni = int(2*getattr(fc, "spin", fractions.Fraction(1, 2)))
        info["E_free"] = lowest(lambda s: s % 2 == ni % 2)
        info["E_screened"] = lowest(lambda s: s % 2 == (ni + 1) % 2)
        if info["E_free"] is not None and info["E_screened"] is not None:
            info["E_ysr"] = info["E_free"] - info["E_screened"]
    if not dmrg: # ED: the plain ground state has a definite parity
        return copy.deepcopy(fc), info
    c = chains.pop(szg)
    chains.clear()
    c.promote_to_dense()
    return c, info
