"""Shared check for the julia_live site types (mpsjulialive/sites.jl and
sites_custom.jl): ground-state energy and a few vev() against mode="ED",
over more than one random-MPS start."""
import pytest

E_TOL = 1e-6  # DMRG energy against ED
V_TOL = 1e-5  # DMRG vev and correlators against ED


def check_julia_live_against_ed(chain, h, observables, maxm=40, nsweeps=20,
        starts=2):
    """Solve `h` on julia_live `starts` times, each a cold random start
    (set_hamiltonian restarts the chain), and compare the energy and
    every observable against ED. Returns the ED energy."""
    chain.set_hamiltonian(h)
    e_ed = chain.gs_energy(mode="ED")
    v_ed = [chain.vev(o, mode="ED") for o in observables]
    chain.setup_julia()
    chain.maxm = maxm
    chain.nsweeps = nsweeps
    for start in range(starts):
        chain.set_hamiltonian(h) # cold start
        e = chain.gs_energy(mode="DMRG")
        assert e == pytest.approx(e_ed, abs=E_TOL), start
        for k, o in enumerate(observables):
            v = chain.vev(o, mode="DMRG")
            assert abs(v - v_ed[k]) < V_TOL, (start, k, v, v_ed[k])
    return e_ed
