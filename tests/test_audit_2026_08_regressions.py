"""Regression tests for the bugs found by the 2026-08 cross-backend audit.

Every test here locks in one finding from `docs/audit_2026_08_hole_hunt.md`,
which records the original symptom, the reproduction and the reviewer's
analysis for each. They are gathered in one file (rather than spread over
the topical ones) because they share a shape: almost all of them are about
a call *silently* doing something other than what was asked, so each test
asserts that the wrong answer is no longer returned, or that what used to
be silent now raises.

Chains are deliberately tiny -- ED is exact at this size and is the
reference wherever one is needed.
"""

import numpy as np
import pytest

from dmrgpy import spinchain, fermionchain, timedependent, thermal


# ---------------------------------------------------------------- helpers

def heisenberg(n=6, itensor_version=3, field=0.0, maxm=30, nsweeps=10):
    """Uniform S=1/2 Heisenberg chain, optionally with an Sz field."""
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=itensor_version)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
              + sc.Sz[i] * sc.Sz[i + 1]
    if field: h = h + field * sc.Sz[0]
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps = maxm, nsweeps
    return sc


def hopping_chain(n=4, itensor_version=3):
    """Spinless nearest-neighbour hopping chain."""
    fc = fermionchain.Fermionic_Chain(n, itensor_version=itensor_version)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
    fc.set_hamiltonian(h)
    fc.maxm, fc.nsweeps = 20, 8
    return fc


# --------------------------------------------- solver parameters are read

@pytest.mark.parametrize("itensor_version", [3, "python"])
def test_gs_energy_reruns_when_maxm_changes(itensor_version):
    """The textbook convergence ramp on a single chain.

    gs_energy()/get_gs() short-circuited on `computed_gs` alone, i.e.
    before groundstate.gs_energy_single() and its send-cache (which does
    key on maxm/nsweeps/cutoff/noise/ramp/sector) were ever reached, so
    every later call returned the first, least converged energy: measured
    -4.047 at every maxm from 2 to 40 on this chain, against an exact
    -4.258.
    """
    sc = heisenberg(n=6, itensor_version=itensor_version, nsweeps=8)
    exact = heisenberg(n=6, itensor_version=itensor_version).gs_energy(mode="ED")
    es = []
    for m in [2, 4, 10, 20, 40]:
        sc.maxm = m
        es.append(sc.gs_energy())
    assert es[0] > es[-1] + 1e-3, "the ramp never re-solved: %s" % es
    assert es[-1] == pytest.approx(exact, abs=1e-6)
    # the ground state itself has to follow the energy
    wf = sc.get_gs()
    assert wf.dot(sc.hamiltonian * wf).real == pytest.approx(exact, abs=1e-6)


@pytest.mark.parametrize("itensor_version", [3, "python"])
def test_repeated_gs_energy_with_unchanged_parameters_is_cached(itensor_version):
    """...and the fix must not make an unchanged repeat call re-solve."""
    sc = heisenberg(n=6, itensor_version=itensor_version)
    e1 = sc.gs_energy()
    import time
    t0 = time.time(); e2 = sc.gs_energy(); dt = time.time() - t0
    assert e1 == e2
    assert dt < 0.05, "a repeat call with unchanged parameters re-solved"


def test_maxm_must_be_positive():
    """maxm<=0 was validated nowhere: itensor_version=2/3 died with an
    uncatchable SIGABRT inside ITensor, "python" silently returned a
    bond-dimension-1 energy."""
    sc = heisenberg(n=4)
    for bad in [0, -5]:
        with pytest.raises(ValueError):
            sc.maxm = bad
    with pytest.raises(ValueError):
        sc.maxm = "not a number"
    sc.maxm = 20  # still accepts a good one


# ------------------------------------------------- time evolution methods

def test_unrecognized_tevol_method_raises():
    """A typo'd tevol_method used to run the legacy MPO-Taylor integrator
    silently -- measured 3.5e-2 error against ED where the requested TDVP
    gives 5.9e-5, bit-identical to an explicit "MPO"."""
    sc = heisenberg(n=4)
    for bad in ["tdvp", "TVDP", "TEBD ", "typo", None]:
        sc.tevol_method = bad
        with pytest.raises(ValueError):
            timedependent.evolve_and_measure_dmrg(
                sc, operator=sc.Sz[0], nt=2, dt=0.1, wf=sc.get_gs())
    for good in ["TDVP", "TEBD", "MPO", "AUTO", "TDVP_GSE"]:
        sc.tevol_method = good  # recognized: no exception from the check


@pytest.mark.parametrize("tevol_method", ["TDVP", "TDVP_GSE"])
def test_v3_evolve_and_measure_preserves_the_input_norm(tevol_method):
    """mpscpp3's evolve_and_measure_tdvp/_gse called a normalizing
    tdvp_step without restoring the input norm, so evolution_ABA -- whose
    start state is A|gs>, norm 1/4 here -- came back scaled by 1/||wf||^2
    for every t>0 while C(0) stayed right."""
    sc = heisenberg(n=5, field=0.4, maxm=20, nsweeps=14)
    ref = heisenberg(n=5, field=0.4, maxm=20, nsweeps=14)
    _t, ced = timedependent.evolution_ABA(ref, A=ref.Sx[0], B=ref.Sz[2],
                                          mode="ED", nt=6, dt=0.2)
    sc.tevol_method = tevol_method
    _t, cs = timedependent.evolution_ABA(sc, A=sc.Sx[0], B=sc.Sz[2],
                                         nt=6, dt=0.2)
    ced, cs = np.asarray(ced).real, np.asarray(cs).real
    assert np.max(np.abs(cs - ced)) < 1e-5


# ------------------------------------------------- dynamical correlators

def test_td_submode_satisfies_the_spectral_sum_rule():
    """_fourier_transform_correlator omitted the codebase-wide 1/pi, so
    submode="TD" (and "TDZ", which shares the same tail) returned
    pi*S(w) + Re C(0)*dt/2: the exact sum rule int S(w) dw = <A B> came
    out at 1.28 against 0.25. Uniform in w, so no peak position or width
    could see it."""
    sc = heisenberg(n=4, field=0.3, maxm=20, nsweeps=12)
    name = (sc.Sz[0], sc.Sz[0])
    es = np.linspace(-20., 20., 4001)
    exact = sc.vev(sc.Sz[0] * sc.Sz[0], mode="ED").real  # = 0.25
    integrals = {}
    for submode in ["KPM", "CVM", "TD"]:
        x, y = sc.get_dynamical_correlator(name=name, submode=submode,
                                           es=es, delta=0.1)
        integrals[submode] = np.trapezoid(np.asarray(y).real, np.asarray(x).real)
    assert integrals["KPM"] == pytest.approx(exact, abs=0.02)
    assert integrals["TD"] == pytest.approx(exact, abs=0.05), \
        "TD weight %s against exact %s" % (integrals["TD"], exact)


def test_kpm_rejects_unknown_keywords():
    """submode="KPM" absorbed every unknown keyword into a bare **kwargs
    and discarded it -- including `deconvolve` and `n`, which sat in its
    own signature and did nothing. CVM/CVM_explicit always raised."""
    sc = heisenberg(n=4)
    name = (sc.Sz[0], sc.Sz[0])
    es = np.linspace(0, 3, 40)
    for bad in [{"deconvolve": "pm"}, {"totally_bogus": 42}, {"nsweeps": 1}]:
        with pytest.raises(TypeError):
            sc.get_dynamical_correlator(name=name, es=es, **bad)


def test_unrecognized_submode_raises_valueerror():
    """The dispatch ended in a bare `raise`, surfacing as the unhelpful
    "RuntimeError: No active exception to reraise"."""
    sc = heisenberg(n=4)
    with pytest.raises(ValueError):
        sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]),
                                     submode="NOT_A_SUBMODE",
                                     es=np.linspace(0, 3, 20))


@pytest.mark.parametrize("mode", ["DMRG", "ED"])
def test_string_operator_names_work_on_every_submode(mode):
    """The documented string form (name="ZZ" with i=/j=) worked only for
    the TD/CVM/TDZ submodes under mode="DMRG"; the default submode="KPM",
    submode="EX" and every mode="ED" path died on it several frames deep
    (KPM with a bare `raise`, i.e. "No active exception to reraise")."""
    sc = heisenberg(n=4)
    es = np.linspace(0, 3, 60)
    x, y = sc.get_dynamical_correlator(mode=mode, name="ZZ", i=0, j=0,
                                        es=es, delta=0.2)
    assert np.max(np.abs(np.asarray(y))) > 0.


def test_missing_optional_submodes_raise_instead_of_exiting():
    """submode="CVMimag"/"maxent" import optional modules that are not
    part of this package; the import failure used to print "Not functional
    yet" and call exit(), terminating the caller's process from inside a
    library."""
    sc = heisenberg(n=4)
    es = np.linspace(0, 3, 20)
    for submode in ["CVMimag", "maxent"]:
        with pytest.raises(NotImplementedError):
            sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]),
                                         submode=submode, es=es)


def test_kpm_energy_truncation_is_not_silently_ignored_on_v2():
    """itensor_version=2 implements no energy truncation at all; the flag
    was accepted and did nothing, so the run looked truncated and was
    not."""
    sc = heisenberg(n=4, itensor_version=2)
    sc.kpm_energy_truncate = True
    with pytest.raises(NotImplementedError):
        sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]),
                                     es=np.linspace(0, 3, 20))


# --------------------------------------------- non-Hermitian Hamiltonians

def nh_chain(n=4, itensor_version=3):
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version=itensor_version)
    h = 0
    for i in range(n - 1):
        h = h + sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1] \
              + sc.Sz[i] * sc.Sz[i + 1]
    sc.set_hamiltonian(h + 0.3j * sc.Sz[0])
    sc.maxm, sc.nsweeps = 20, 6
    return sc


def test_non_hermitian_submode_is_no_longer_discarded():
    """The Hermiticity check ran before the submode dispatch and returned
    the explicit CVM resolvent for everything except "KPM", so a
    non-Hermitian Hamiltonian made submode= a no-op: EX/maxent/ROOTN/TD/
    CVM all returned the same curve. EX is documented (ROADMAP.md) as
    non-Hermitian-capable and now actually runs."""
    es = np.linspace(0.2, 3., 12)
    sc = nh_chain()
    _x, y_ex = sc.get_dynamical_correlator(name=[sc.Sz[0], sc.Sz[0]],
                                            submode="EX", es=es, delta=0.1,
                                            nex=4)
    sc2 = nh_chain()
    _x, y_cvm = sc2.get_dynamical_correlator(name=[sc2.Sz[0], sc2.Sz[0]],
                                              submode="CVM_explicit",
                                              es=es, delta=0.1)
    assert not np.allclose(np.asarray(y_ex), np.asarray(y_cvm), atol=1e-8), \
        "EX still returns the CVM_explicit result"


@pytest.mark.parametrize("submode", ["TD", "ROOTN"])
def test_hermitian_only_submodes_raise_on_a_non_hermitian_hamiltonian(submode):
    """...and the ones with no non-Hermitian implementation say so rather
    than quietly answering with a different algorithm."""
    sc = nh_chain()
    with pytest.raises(NotImplementedError):
        sc.get_dynamical_correlator(name=[sc.Sz[0], sc.Sz[0]],
                                     submode=submode,
                                     es=np.linspace(0.2, 3., 12), delta=0.1)


def test_ed_lehmann_reference_is_not_substituted():
    """The worst case of the same bug: mode="ED", submode="ED" is the
    exact reference a user cross-validates against, and it silently became
    the approximate resolvent -- agreeing with it bit for bit."""
    sc = nh_chain()
    with pytest.raises(NotImplementedError):
        sc.get_dynamical_correlator(mode="ED", name=[sc.Sz[0], sc.Sz[0]],
                                     submode="ED",
                                     es=np.linspace(0.2, 3., 12), delta=0.1)


def test_excited_states_beyond_the_hilbert_space_raise():
    """Asking for more states than exist has no answer; the non-Hermitian
    Krylov route spun for minutes printing "state is not normalizable"
    instead of saying so (dcex's own default nex=20 on a 16-dimensional
    chain)."""
    sc = heisenberg(n=4)
    with pytest.raises(ValueError):
        sc.get_excited_states(n=20)  # dimension is 2**4 = 16


# -------------------------------------------------- conserved-sector mode

def test_sector_reduced_dm_is_a_density_matrix():
    """mpscpp3's reduced_dm flattened rho with rho.visit(), which under
    QN-carrying indices enumerates only the stored blocks; bindings.cc
    copied that short vector into an uninitialized (dim,dim) array, so
    get_rdm() returned uninitialized heap memory -- non-Hermitian, and
    different from one run to the next."""
    for sites, field in [(["S=1/2"] * 6, 0.0), (["S=1"] * 4, 0.3)]:
        n = len(sites)
        def build():
            s = spinchain.Spin_Chain(sites)
            h = 0
            for i in range(n - 1):
                h = h + s.Sx[i] * s.Sx[i + 1] + s.Sy[i] * s.Sy[i + 1] \
                      + s.Sz[i] * s.Sz[i + 1]
            if field: h = h + field * s.Sz[0]
            s.set_hamiltonian(h)
            s.maxm, s.nsweeps = 30, 10
            return s
        dense = build().get_rdm(i=1)
        s2 = build()
        s2.set_conserved_sector(Sz=0)
        s2.gs_energy()
        sec = s2.get_rdm(i=1)
        assert np.trace(sec).real == pytest.approx(1.0, abs=1e-6)
        assert np.allclose(sec, sec.conj().T, atol=1e-8)
        assert np.allclose(np.diag(sec).real, np.diag(dense).real, atol=1e-6)


def test_sector_same_site_operator_product_does_not_abort():
    """sector_terms validated each term's NET charge only, so a
    net-neutral term piling an impossible charge on ONE site (Cdag0 C1
    Cdag0 C1, produced routinely by the four-point correlator) passed the
    guard and aborted the whole process inside AutoMPO. Such a product is
    zero by Pauli exclusion, which is what ED answers."""
    fc = hopping_chain(n=4)
    fc.set_conserved_sector(Nf=2)
    fc.gs_energy()
    op = fc.Cdag[0] * fc.C[1] * fc.Cdag[0] * fc.C[1]
    assert abs(complex(fc.vev(op))) < 1e-10
    # ordinary repeated-site products must keep working
    for a, b, c, d in [(0, 1, 1, 0), (1, 0, 0, 1), (0, 0, 1, 1)]:
        o = fc.Cdag[a] * fc.C[b] * fc.Cdag[c] * fc.C[d]
        assert np.isfinite(abs(complex(fc.vev(o))))


@pytest.mark.parametrize("itensor_version", [3, "python"])
def test_sector_site_and_pair_entropies_work(itensor_version):
    """reduced_dm_projective applied a charge-raising projector (S+) to
    the state, which conserved-sector mode rejects, so site/pair entropy
    and mutual information all raised there. Evaluated as a single
    <wf|Pa^dag Pb|wf> sandwich they conserve the charge (and the
    cross-charge entries are exactly zero in a fixed-charge state)."""
    ref = heisenberg(n=6, itensor_version=itensor_version)
    wf = ref.get_gs()
    s_ref = ref.get_site_entropy(wf, 2)
    mi_ref = ref.get_mutual_information(wf, 1, 4)
    sec = heisenberg(n=6, itensor_version=itensor_version)
    sec.set_conserved_sector(Sz=0)
    wf2 = sec.get_gs()
    assert sec.get_site_entropy(wf2, 2) == pytest.approx(s_ref, abs=1e-5)
    assert sec.get_mutual_information(wf2, 1, 4) == pytest.approx(mi_ref, abs=1e-5)
    sec.get_pair_entropy(wf2, 1, 4)  # must not raise


def test_sector_correlation_matrix_default_mode_works():
    """get_correlation_matrix's hardcoded default dmmode="fast" applies
    single-fermion operators to the state, which changes the particle
    number and so always raised in sector mode -- while
    dmmode="explicit"/"full" returned the right answer on the same
    chain."""
    fc = hopping_chain(n=4)
    ref = fc.get_correlation_matrix()
    fc2 = hopping_chain(n=4)
    fc2.set_conserved_sector(Nf=2)
    fc2.gs_energy()
    got = fc2.get_correlation_matrix()  # default dmmode
    assert np.allclose(np.asarray(got), np.asarray(ref), atol=1e-5)


# ------------------------------------------------------ pyitensor backend

def test_pyitensor_rejects_an_out_of_range_site():
    """An operator on a site outside the chain fell through pyitensor's
    per-site loop and silently became the identity: vev(Sz[10]) on a
    6-site chain returned 1.0, and 5.0*Sz[10] in the Hamiltonian shifted
    the energy by exactly +5. v2, v3 and ED all raise IndexError."""
    sc = heisenberg(n=6, itensor_version="python")
    with pytest.raises(IndexError):
        sc.vev(sc.get_operator("Sz", 10))
    with pytest.raises(IndexError):
        sc.vev(sc.get_operator("Sz", -1))


def test_pyitensor_reduced_dm_at_the_last_site():
    """pyitensor's reduced_dm indexed psi.A(site+1) unconditionally, so
    the last site raised IndexError while v2/v3 returned the matrix."""
    sc = heisenberg(n=4, itensor_version="python")
    ref = heisenberg(n=4, itensor_version=3)
    for site in [0, 3]:
        got = np.diag(sc.get_rdm(i=site)).real
        want = np.diag(ref.get_rdm(i=site)).real
        assert np.allclose(np.sort(got), np.sort(want), atol=1e-5)


# ------------------------------------------------------------ odds & ends

def test_same_site_fermion_product_is_zero_everywhere():
    """C_i C_i = 0, but jordanwigner.CC(i,i) returned the Python int 0, so
    every DMRG backend raised "'float' object has no attribute 'op'" for a
    legitimate vev(C[i]*C[i]) that ED answered with 0j."""
    fc = hopping_chain(n=4)
    for mode in ["DMRG", "ED"]:
        assert abs(complex(fc.vev(fc.C[0] * fc.C[0], mode=mode))) < 1e-12
        assert abs(complex(fc.vev(fc.Cdag[0] * fc.Cdag[0], mode=mode))) < 1e-12


def test_thermal_chain_honors_its_constructor_kwargs():
    """Thermal_Spin_Chain accepted **kwargs and dropped them, so
    itensor_version= silently built the default backend instead."""
    tc = thermal.Thermal_Spin_Chain(["S=1/2"] * 2, T=0.5,
                                     itensor_version="python")
    assert tc.MBChain.itensor_version == "python"


def test_thermal_chain_rejects_a_negative_temperature():
    """T<0 was silently treated as T=0 (and T exactly 1e-5 fell through
    both branches into an UnboundLocalError)."""
    with pytest.raises(ValueError):
        thermal.Thermal_Spin_Chain(["S=1/2"] * 2, T=-1.0)


def test_get_distribution_dispatches_through_get_mode():
    """These were the only public entry points that never called
    get_mode(), so every route into ED handed session-only code an ED
    State: "'State' object has no attribute 'cpp_handle'"."""
    sc = heisenberg(n=4)
    with pytest.raises(NotImplementedError):
        sc.get_distribution_moments(mode="ED", X=sc.Sz[0])


def test_dead_evolution_method_is_gone():
    """Many_Body_Chain.evolution() called timedependent.evolution(), which
    does not exist, so it raised AttributeError on every backend."""
    sc = heisenberg(n=4)
    assert not hasattr(sc, "evolution")
