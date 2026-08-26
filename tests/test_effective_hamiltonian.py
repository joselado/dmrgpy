"""Coverage for effectivehamiltonian.py.

The reference case throughout is the Hubbard dimer at strong coupling,
which maps onto a two-spin exchange model with J = 4t^2/U. Adding a
Peierls phase phi twists the two spin species oppositely, so the exchange
keeps its magnitude while its XY part rotates by 2*phi:

    H_eff = J [ cos(2 phi)(SxSx + SySy) + sin(2 phi)(SxSy - SySx) + SzSz ]

That gives an analytic target to check the fit against, rather than a
golden value copied out of a previous run.

Several of these are regression tests for defects the module had:
`n` cutting a degenerate multiplet returned meaningless couplings with no
warning; the fit used a random-start nonlinear optimizer on an
exactly-linear problem, so it was neither accurate nor reproducible;
`get_projection_operators` hardcoded self.ns//2 and so covered only half
the sites of any non-spinful-fermionic chain; Spin_Chain's own entry
point raised TypeError; and `operators=`/`tol=` were declared but ignored.
"""
import numpy as np
import pytest

from dmrgpy import fermionchain, spinchain
from dmrgpy import effectivehamiltonian as eh

# h = h + h.get_dagger() below doubles the (already Hermitian) Hubbard
# term while leaving |t| = 1, so the effective U is 40, not 20
U_EFF = 40.0
J_EXACT = 4.0*1.0**2/U_EFF


def hubbard_dimer(phi=0.0):
    """Two spinful fermionic sites, strong Hubbard U, Peierls phase phi
    (in units of pi) on the hopping."""
    n = 2
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    h = 0
    t = 1.0*np.exp(1j*phi*np.pi)
    for i in range(n-1):
        h = h + t*fc.Cdagup[i]*fc.Cup[i+1]
        h = h + np.conjugate(t)*fc.Cdagdn[i]*fc.Cdn[i+1]
    for i in range(n):
        h = h + 0.5*U_EFF*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
    h = h + h.get_dagger()
    fc.set_hamiltonian(h)
    return fc


def heisenberg_spin_chain(n=4):
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1]
        h = h + sc.Sy[i]*sc.Sy[i+1]
        h = h + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    return sc


def test_hubbard_dimer_gives_isotropic_heisenberg_exchange():
    """With no Peierls phase the effective model is plain Heisenberg:
    equal xx/yy/zz couplings of 4t^2/U, and nothing else."""
    coef = eh.get_effective_hamiltonian_couplings(hubbard_dimer(),
            method="single", mode="ED")
    for pair in [("S^x_1", "S^x_2"), ("S^y_1", "S^y_2"), ("S^z_1", "S^z_2")]:
        assert np.real(coef[pair]) == pytest.approx(J_EXACT, abs=1e-3)
    # no spurious XY rotation at phi=0
    assert abs(np.real(coef.get(("S^x_1", "S^y_2"), 0.0))) < 1e-6


@pytest.mark.parametrize("phi", [0.0, 0.1, 0.2, 0.3, 0.5])
def test_peierls_phase_rotates_the_xy_exchange_by_twice_the_phase(phi):
    """The analytic content of this module's own example: the fitted
    couplings must follow J*cos(2 phi) / J*sin(2 phi) / J."""
    coef = eh.get_effective_hamiltonian_couplings(hubbard_dimer(phi),
            method="single", mode="ED")
    xx = np.real(coef.get(("S^x_1", "S^x_2"), 0.0))
    xy = np.real(coef.get(("S^x_1", "S^y_2"), 0.0))
    yx = np.real(coef.get(("S^y_1", "S^x_2"), 0.0))
    zz = np.real(coef.get(("S^z_1", "S^z_2"), 0.0))

    assert xx == pytest.approx(J_EXACT*np.cos(2*np.pi*phi), abs=1e-3)
    assert xy == pytest.approx(J_EXACT*np.sin(2*np.pi*phi), abs=1e-3)
    assert yx == pytest.approx(-J_EXACT*np.sin(2*np.pi*phi), abs=1e-3)
    assert zz == pytest.approx(J_EXACT, abs=1e-3)
    # the phase only rotates the exchange, it does not change its size
    assert np.sqrt(xx**2 + xy**2) == pytest.approx(J_EXACT, abs=1e-3)


def test_full_and_single_methods_agree_on_this_manifold():
    """P A P B P and P (A B) P coincide when the manifold is closed under
    the operators, which it is for the dimer's singlet+triplet."""
    a = eh.get_effective_hamiltonian_couplings(hubbard_dimer(0.1),
            method="single", mode="ED")
    b = eh.get_effective_hamiltonian_couplings(hubbard_dimer(0.1),
            method="full", mode="ED")
    assert set(a) == set(b)
    for key in a:
        assert np.real(a[key]) == pytest.approx(np.real(b[key]), abs=1e-3)


def test_unknown_method_is_rejected():
    with pytest.raises(ValueError):
        eh.get_effective_hamiltonian_couplings(hubbard_dimer(),
                method="nonsense", mode="ED")


@pytest.mark.parametrize("n", [2, 3])
def test_cutting_a_degenerate_multiplet_raises(n):
    """The dimer's low manifold is a singlet plus a 3-fold degenerate
    triplet, so only n=1 and n=4 cut cleanly. n=2 or n=3 leaves a
    manifold that is not symmetry-invariant; the fit used to return
    plausible-looking couplings ~50x the correct ones with no warning.
    """
    with pytest.raises(ValueError, match="degenerate multiplet"):
        eh.get_effective_hamiltonian_couplings(hubbard_dimer(),
                method="single", mode="ED", n=n)


def test_degeneracy_error_names_a_usable_n():
    """The error has to be actionable: it reports which n do cut
    cleanly among the states it looked at."""
    with pytest.raises(ValueError) as excinfo:
        eh.get_effective_hamiltonian_couplings(hubbard_dimer(),
                method="single", mode="ED", n=3)
    assert "cut cleanly" in str(excinfo.value)


def test_fit_is_reproducible():
    """fit_matrix used to minimize from an unseeded random start, so two
    identical calls disagreed at ~2.6e-5 -- enough for a coupling near
    the cutoff to be returned by one call and dropped by the next."""
    a = eh.get_effective_hamiltonian_couplings(hubbard_dimer(0.1),
            method="single", mode="ED")
    b = eh.get_effective_hamiltonian_couplings(hubbard_dimer(0.1),
            method="single", mode="ED")
    assert set(a) == set(b)
    for key in a:
        assert a[key] == b[key]  # bit for bit, not merely close


def test_fit_matrix_solves_the_least_squares_problem_exactly():
    """A linear problem solved linearly: the residual must be at machine
    precision, not at an optimizer's convergence tolerance (the old
    random-start minimize left 3.5e-6 behind)."""
    rng = np.random.default_rng(0)
    basis = {}
    for k in range(4):
        m = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4))
        basis["op%d" % k] = m
    x_true = np.array([0.3, -1.2, 2.5, 0.75])
    h = sum(x_true[k]*basis["op%d" % k] for k in range(4))

    coef = eh.fit_matrix(h, basis, cutoff=1e-12)
    got = np.array([coef["op%d" % k] for k in range(4)])
    assert got == pytest.approx(x_true, abs=1e-10)


def test_fit_matrix_raises_on_a_linearly_dependent_basis():
    """A degenerate basis has no unique answer -- the coefficients can be
    shifted by any null-space vector without changing the residual -- so
    the fit refuses rather than reporting one arbitrary representative."""
    rng = np.random.default_rng(1)
    a = rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3))
    b = rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3))
    basis = {"a": a, "b": b, "a_plus_b": a + b}  # rank 2, three operators
    with pytest.raises(ValueError, match="linearly dependent"):
        eh.fit_matrix(a, basis)


def test_projection_operators_cover_every_site_of_a_spin_chain():
    """get_projection_operators hardcoded self.ns//2, which is right only
    for Spinful_Fermionic_Chain (two spinless sites per orbital): on a
    4-site Spin_Chain it silently returned operators for sites 1-2 only,
    so the fit ran against a basis that could not represent H."""
    sc = heisenberg_spin_chain(n=4)
    op = eh.get_projection_operators(sc)
    assert len(op) == 3*sc.ns
    assert sorted({k.split("_")[1] for k in op}) == ["1", "2", "3", "4"]

    fc = hubbard_dimer()
    # ...and one entry per *orbital* for a spinful fermionic chain, whose
    # self.ns counts the interleaved up/down sites
    assert len(eh.get_projection_operators(fc)) == 3*(fc.ns//2)


def test_spin_chain_entry_point_runs():
    """Spin_Chain.get_effective_hamiltonian() used to raise TypeError
    unconditionally: spinchain.py passed a name="XX" that nothing
    downstream accepted."""
    sc = heisenberg_spin_chain(n=4)
    # n=1 is the ground state alone, which always cuts cleanly here
    out = sc.get_effective_hamiltonian(mode="ED", n=1)
    assert isinstance(out, str)


def test_get_heff_recovers_the_exchange_couplings():
    """The operators= path: fitting with the three exchange operators
    directly must give three couplings equal to 4t^2/U."""
    fc = hubbard_dimer()
    ops = [fc.Sx[0]*fc.Sx[1], fc.Sy[0]*fc.Sy[1], fc.Sz[0]*fc.Sz[1]]
    coef = np.real(fc.get_heff(mode="ED", operators=ops))
    assert len(coef) == len(ops)
    assert coef == pytest.approx([J_EXACT]*3, abs=1e-3)


def test_get_heff_returns_one_coefficient_per_operator():
    """A coefficient below the cutoff must come back as 0.0 rather than
    shortening the list and silently misaligning it against `operators`."""
    fc = hubbard_dimer()
    ops = [fc.Sx[0]*fc.Sx[1], fc.Sy[0]*fc.Sy[1], fc.Sz[0]*fc.Sz[1],
           fc.Sx[0]]  # the last one carries no weight in H_eff
    coef = np.real(fc.get_heff(mode="ED", operators=ops))
    assert len(coef) == 4
    assert coef[3] == 0.0


def test_get_heff_without_operators_raises():
    """It used to return the NotImplemented singleton, which fails
    obscurely somewhere downstream instead of here."""
    with pytest.raises(ValueError):
        hubbard_dimer().get_heff(mode="ED")


def test_operators_argument_is_honoured_by_couplings():
    """`operators=` was declared and never read: the function always used
    get_projection_operators, so a caller-supplied basis was silently
    replaced."""
    fc = hubbard_dimer()
    op = {"S^z_1": fc.Sz[0], "S^z_2": fc.Sz[1]}
    coef = eh.get_effective_hamiltonian_couplings(fc, method="single",
            mode="ED", operators=op)
    # only the given operators can appear in the keys
    for key in coef:
        for name in key:
            assert name in op
    assert np.real(coef[("S^z_1", "S^z_2")]) == pytest.approx(J_EXACT, abs=1e-3)


def test_tol_argument_is_honoured():
    """`tol=` was declared and never read. Raising it above a coupling's
    magnitude must drop that coupling."""
    fc = hubbard_dimer()
    fine = eh.get_effective_hamiltonian_couplings(fc, method="single",
            mode="ED", tol=1e-4)
    coarse = eh.get_effective_hamiltonian_couplings(fc, method="single",
            mode="ED", tol=1.0)  # far above every exchange coupling
    assert len(fine) > 0
    assert len(coarse) < len(fine)


def test_dmrg_and_ed_agree():
    """The DMRG path goes through excited states that are only
    approximately orthogonal, so it exercises the Lowdin orthonormalization
    in a way the ED path does not.

    nsweeps/maxm are raised above the chain defaults because this test was
    intermittently failing: the DMRG side starts from a random MPS, so its
    result depends on the global RNG state, which depends on which tests ran
    before it -- it passed every time in isolation and failed roughly one
    run in forty inside the full suite. Measured worst-case |ED-DMRG| over
    many seeds, which shows it is under-convergence rather than a tolerance
    that is merely unlucky:

        chain defaults          8.5e-05   (vs the 1e-4 assertion below)
        nsweeps=20              4.6e-05
        nsweeps=40              3.3e-06
        nsweeps=40, maxm=120    2.3e-07

    So the fix is to converge the excited states rather than to loosen the
    assertion (an earlier attempt raised the tolerance to 5e-4; a 60-seed
    scan then found a 6.6e-4 outlier, i.e. tolerance-chasing was losing to
    the tail). The dimer is two sites, so 40 sweeps costs almost nothing.
    """
    fc_ed, fc_dmrg = hubbard_dimer(0.1), hubbard_dimer(0.1)
    for fc in (fc_ed, fc_dmrg):
        fc.nsweeps = 40
        fc.maxm = 120
    a = eh.get_effective_hamiltonian_couplings(fc_ed, method="single", mode="ED")
    b = eh.get_effective_hamiltonian_couplings(fc_dmrg, method="single", mode="DMRG")
    assert set(a) == set(b)
    for key in a:
        assert np.real(a[key]) == pytest.approx(np.real(b[key]), abs=1e-4)


def test_orthonormalize_fixes_a_non_orthonormal_basis():
    """Directly: Lowdin orthonormalization of a deliberately skewed basis
    must return a basis with an identity Gram matrix that spans the same
    space."""
    fc = hubbard_dimer()
    _, ws = fc.get_excited_states(n=4, mode="ED")
    ws = list(ws)
    skewed = [ws[0], ws[1] + 0.3*ws[0], ws[2] - 0.2*ws[1], ws[3]]
    S = np.array([[skewed[i].overlap(skewed[j]) for j in range(4)]
                  for i in range(4)])
    assert np.max(np.abs(S - np.identity(4))) > 0.1  # genuinely skewed

    fixed = eh.orthonormalize(skewed)
    S2 = np.array([[fixed[i].overlap(fixed[j]) for j in range(4)]
                   for i in range(4)])
    assert np.max(np.abs(S2 - np.identity(4))) < 1e-10


def test_dict2latex_formats_signs_and_terminators():
    """The output used to put 'H =' and the prefactor outside the math
    block, and to terminate every term with '+' -- including the last,
    leaving a dangling '+' before the closing bracket."""
    out = eh.dict2latex({("S^x_1", "S^x_2"): 1.0, ("S^y_1", "S^x_2"): -0.5})
    assert out.startswith("\\[")
    assert out.rstrip().endswith("\\]")
    assert "+  -" not in out       # negatives use a '-' separator
    assert "-  0.5" in out
    assert "+ \n\\right)" not in out  # no dangling terminator


def test_dict2latex_handles_an_empty_or_zero_fit():
    """np.max over an empty list, and a division by a zero prefactor,
    were both reachable."""
    assert "0" in eh.dict2latex({})
    assert "0" in eh.dict2latex({("S^x_1", "S^x_2"): 0.0})
