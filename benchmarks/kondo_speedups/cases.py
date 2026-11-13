"""The calls that docs/kondo_speedup_plan.md times. Each case has a setup
(the chain, and on DMRG its solved ground state, which is not timed), the
public call, the items of the plan it measures, and the functions whose
cumulative time the profile reports. The chains are those of the records:
docs/audit_2026_11_hole_hunt.md (ED and TASK),
docs/new_methods/kondo_krylov_thermal/README.md section 6 (Krylov) and
docs/new_methods/kondo_leading_log/README.md section 7 (leading logs)."""
import numpy as np

from dmrgpy import spinchain

ED_WATCH = ["edkondo.py:__init__", "edkondo.py:_eigenbasis_spin_operator",
            "conductance.py:second_order_dIdV",
            "conductance.py:third_order_kondo_dIdV",
            "conductance.py:third_order_potential_dIdV",
            "conductance.py:_triple_product_coefficients",
            "conductance.py:_levi_civita", "conductance.py:_kondo_hub_blocks",
            "edkondo.py:eigenpairs",
            "stepfunctions.py:__call__", "stepfunctions.py:F0", "eigh"]
TASK_WATCH = ["taskkondo.py:__init__", "task.py:projected_spectrum",
              "task.py:block_krylov", "task.py:apply_many",
              "conductance.py:second_order_dIdV",
              "conductance.py:third_order_kondo_dIdV", "get_excited_states"]
KRYLOV_WATCH = ["krylovkondo.py:__init__", "krylovkondo.py:low_states",
                "krylovkondo.py:state_spectrum", "krylovkondo.py:block_krylov",
                "conductance.py:second_order_dIdV",
                "conductance.py:third_order_kondo_dIdV",
                "conductance.py:third_order_potential_dIdV",
                "conductance.py:_levi_civita"]
LL_WATCH = ["edkondo.py:__init__", "leadinglog.py:leading_log_dIdV",
            "leadinglog.py:_joint_flow", "leadinglog.py:tip_vertex",
            "leadinglog.py:_hole_sum_table", "leadinglog.py:_hole_sum",
            "thermalkernel.py:__call__", "thermalkernel.py:__init__"]


def heisenberg_ed(n):
    """The ED chain of the 2026-11 hunt: J=4 meV, fields on sites 0 and 2."""
    sc = spinchain.Spin_Chain(["1/2"]*n)
    h = 0
    for i in range(n - 1):
        h = h + 0.004*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    sc.set_hamiltonian(h + 0.001*sc.Sz[0] + 0.0003*sc.Sx[2])
    return sc


def ed_call(n, T, nb):
    def setup():
        return heisenberg_ed(n)
    def call(sc):
        eV = np.linspace(-0.01, 0.01, nb)
        return sc.get_kondo_spectrum(eV, site=0, Jrho_s=-0.05, U=0.25, T=T,
                                     mode="ED")[1]
    return setup, call


def disordered_xxz(n=20):
    """The 20-site TASK chain of the 2026-11 hunt, ground state solved."""
    rng = np.random.default_rng(5)
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version="python")
    J = 1e-3*rng.uniform(0.5, 1.5, n - 1)
    h = 0
    for i in range(n - 1):
        h = h + J[i]*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + 1.5*sc.Sz[i]*sc.Sz[i+1])
    for i in range(n):
        h = h + 6e-4*(-1)**i*sc.Sz[i]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps, sc.cutoff = 40, 14, 1e-10
    sc.gs_energy()
    return sc


def task_call(order, nkry, T, nb):
    def call(sc):
        eV = np.linspace(-6e-3, 6e-3, nb)
        return sc.get_kondo_spectrum(eV, site=3, Jrho_s=-0.1, U=0.3, T=T,
                                     omega0=20e-3, Gamma0=5e-5, order=order,
                                     mode="DMRG", submode="TASK", nkry=nkry)[1]
    return disordered_xxz, call


def odd_chain():
    """A 9-site Heisenberg chain at maxm=10: a truncated Kramers doublet."""
    sc = spinchain.Spin_Chain(["S=1/2"]*9, itensor_version="python")
    h = 0
    for i in range(8):
        h = h + 1e-3*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = 10, 20
    sc.gs_energy()
    return sc


def ngs_scan(sc):
    """Item 10: a scan of U at n_gs=2, which solves the manifold per call."""
    eV = np.linspace(-4e-3, 4e-3, 41)
    out = [sc.get_kondo_spectrum(eV, site=0, Jrho_s=-0.1, U=U, T=0.,
                                 omega0=20e-3, Gamma0=5e-5, order=3,
                                 mode="DMRG", nkry=60, n_gs=2)[1]
           for U in (0.0, 0.1, 0.2, 0.3)]
    return np.concatenate(out)


def krylov_chain(n=12):
    """Section 6 of the Krylov README: J = 1 meV (1 + 0.3 cos i), a 120
    degree field of 0.2 meV in the xz plane."""
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    h = 0
    for i in range(n - 1):
        J = 1e-3*(1. + 0.3*np.cos(i))
        h = h + J*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    a = 2.*np.pi/3.
    for i in range(n):
        h = h + 0.2e-3*(np.cos(a)*sc.Sz[i] + np.sin(a)*sc.Sx[i])
    sc.set_hamiltonian(h)
    return sc


def krylov_call(sc):
    eV = np.linspace(-4e-3, 4e-3, 61)
    return sc.get_kondo_spectrum(eV, site=0, Jrho_s=-0.05, U=0.25, T=1.0,
                                 mode="ED", submode="KRYLOV", nkry=100)[1]


def ll_chain(n):
    """Section 7 of the leading-log README: J = 1 meV, a field of 0.3 meV
    along z and 0.1 meV along x on every site."""
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    h = 0
    for i in range(n - 1):
        h = h + 1e-3*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    for i in range(n):
        h = h + 0.3e-3*sc.Sz[i] + 0.1e-3*sc.Sx[i]
    sc.set_hamiltonian(h)
    return sc


def ll_call(T, mem_gb=None):
    def call(sc):
        eV = np.linspace(-3e-3, 3e-3, 61) + 1e-5
        kw = {} if mem_gb is None else dict(mem_gb=mem_gb)
        return sc.get_kondo_spectrum(eV, site=0, Jrho_s=-0.05, U=0.25, T=T,
                                     order="LL", **kw)[1]
    return call


def _case(what, setup, call, watch, reps):
    return dict(what=what, setup=setup, call=call, watch=watch, reps=reps)


CASES = {}
for name, (n, T, nb, items) in {
        "ed_T0_n10": (10, 0., 41, "1 2 3 5"),
        "ed_T1K_n10": (10, 1., 41, "1 3"),
        "ed_T20K_n9": (9, 20., 41, "2 3 4 5"),
        "ed_T10K_n10_b241": (10, 10., 241, "4 5"),
        }.items():
    setup, call = ed_call(n, T, nb)
    CASES[name] = _case("mode=ED, %d spins, T=%g K, %d biases; items %s" % (n, T, nb, items),
                        setup, call, ED_WATCH, 3)
for name, (order, nkry, T, nb, items) in {
        "task_order2_nkry100": (2, 100, 0., 121, "6 7 8 9"),
        "task_order3_nkry300_T03": (3, 300, 0.3, 241, "7 8 9"),
        }.items():
    setup, call = task_call(order, nkry, T, nb)
    CASES[name] = _case("TASK on python, 20 sites, order %d, nkry=%d, T=%g K, %d biases; items %s"
                        % (order, nkry, T, nb, items), setup, call, TASK_WATCH, 3)
CASES["task_ngs2_scan"] = _case("TASK on python, 9-site doublet, n_gs=2, 4 values of U; item 10",
                                odd_chain, ngs_scan, TASK_WATCH, 2)
CASES["krylov_n12_T1K"] = _case("submode=KRYLOV, 12 spins, T=1 K, 61 biases; items 2 11 12",
                                krylov_chain, krylov_call, KRYLOV_WATCH, 2)
CASES["ll_n9_T0"] = _case("order=LL on the CPU, 9 spins, T=0, 61 biases; items 13 15",
                          lambda: ll_chain(9), ll_call(0.), LL_WATCH, 2)
CASES["ll_n8_T2K"] = _case("order=LL on the CPU, 8 spins, T=2 K, 61 biases, mem_gb=12; items 13 14",
                           lambda: ll_chain(8), ll_call(2., mem_gb=12.), LL_WATCH, 1)
