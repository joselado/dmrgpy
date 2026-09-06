"""What a chain runs on when nobody names a backend -- in particular on a
machine with no compiled C++ extension, which is every `pip install
dmrgpy`.

The wheel deliberately ships no C++ at all (CLAUDE.md, "Packaging /
PyPI"). Until `cppext.default_backend()` existed, that meant every chain
built with the default backend fell through mode.py's
extension-not-compiled branch to **exact diagonalization** -- so a script
written against the compiled backend kept running and silently changed
algorithm, to one that cannot reach the sizes an MPS solver exists for.
The pure-Python pyitensor backend has no compiled precondition of any
kind, so it is the right default there, and these tests pin that.

They simulate the missing extension by emptying `cppext`'s module-level
backend cache, so they exercise the real decision on a machine where the
extension *is* built.
"""
import pytest

from dmrgpy import cppext, spinchain


@pytest.fixture
def no_cpp_extension(monkeypatch):
    """Make cppext report every compiled backend as unavailable, without
    touching the pure-Python one (which has nothing to be missing)."""
    monkeypatch.setitem(cppext._backends, 2, None)
    monkeypatch.setitem(cppext._backends, 3, None)
    assert not cppext.available(3)
    assert cppext.available("python")


def test_default_backend_is_the_cpp_version_when_it_is_compiled():
    if not cppext.available(cppext.DEFAULT_ITENSOR_VERSION):
        pytest.skip("no compiled C++ extension in this environment")
    assert cppext.default_backend() == cppext.DEFAULT_ITENSOR_VERSION


def test_default_backend_falls_back_to_pyitensor(no_cpp_extension):
    assert cppext.default_backend() == "python"


def test_a_chain_with_no_itensor_version_picks_pyitensor(no_cpp_extension):
    """The behaviour change itself: no kwarg, no C++, and the chain is on a
    real DMRG backend rather than on ED."""
    sc = spinchain.Spin_Chain(["S=1/2"] * 6)
    assert sc.itensor_version == "python"


def test_that_chain_actually_runs_dmrg_and_agrees_with_ed(no_cpp_extension):
    sc = spinchain.Spin_Chain(["S=1/2"] * 6)
    sc.set_hamiltonian(sum(sc.Sx[i] * sc.Sx[i + 1] + sc.Sy[i] * sc.Sy[i + 1]
                           + sc.Sz[i] * sc.Sz[i + 1] for i in range(5)))
    from dmrgpy import mode
    assert mode.get_mode(sc, mode="DMRG") == "DMRG"
    assert sc.gs_energy(mode="DMRG") == pytest.approx(sc.gs_energy(mode="ED"),
                                                      abs=1e-6)


def test_an_explicit_cpp_version_still_falls_back_to_ed(no_cpp_extension):
    """Only the *implicit* choice moved. A caller who named
    itensor_version=3 asked for that backend specifically, so mode.py's
    long-standing ED fallback still applies to them -- silently rerouting
    an explicit request to a different DMRG backend would be a worse
    surprise than the fallback that has always been there."""
    from dmrgpy import mode
    sc = spinchain.Spin_Chain(["S=1/2"] * 6, itensor_version=3)
    assert sc.itensor_version == 3
    assert mode.get_mode(sc, mode="DMRG") == "ED"


def test_explicitly_asking_for_python_is_unchanged(no_cpp_extension):
    sc = spinchain.Spin_Chain(["S=1/2"] * 6, itensor_version="python")
    assert sc.itensor_version == "python"


def test_mode_ed_is_still_a_deliberate_choice(no_cpp_extension):
    """ED is not what the fallback demoted -- it stays exactly what it was,
    the cross-check the rest of tests/ is built on."""
    from dmrgpy import mode
    sc = spinchain.Spin_Chain(["S=1/2"] * 6)
    assert mode.get_mode(sc, mode="ED") == "ED"


@pytest.mark.parametrize("n", [1])
def test_pyitensor_falls_back_to_ed_below_two_sites(n):
    """pyitensor's DMRG is two-site (dmrg.py's `for i in range(1,n)`), so a
    one-site chain has no update to make: the sweep body never runs and
    dmrg() returns the `energy = None` it started with. mode.py routes
    those to ED, the same way it does for ITensor v3 below three sites.
    Before that, gs_energy(mode="DMRG") returned None here."""
    from dmrgpy import mode
    sc = spinchain.Spin_Chain(["S=1/2"] * n, itensor_version="python")
    sc.set_hamiltonian(0.8 * sc.Sz[0] + 0.6 * sc.Sx[0])
    assert mode.get_mode(sc, mode="DMRG") == "ED"
    assert sc.gs_energy(mode="DMRG") == pytest.approx(-0.5, abs=1e-8)


def test_two_sites_is_still_real_dmrg_on_pyitensor():
    """The cutoff is 2, not v3's 3: a two-site chain has exactly one
    two-site update and solves correctly."""
    from dmrgpy import mode
    sc = spinchain.Spin_Chain(["S=1/2"] * 2, itensor_version="python")
    sc.set_hamiltonian(sc.Sx[0] * sc.Sx[1] + sc.Sy[0] * sc.Sy[1]
                       + sc.Sz[0] * sc.Sz[1] + 0.3 * sc.Sz[0])
    assert mode.get_mode(sc, mode="DMRG") == "DMRG"
    assert sc.gs_energy(mode="DMRG") == pytest.approx(sc.gs_energy(mode="ED"),
                                                      abs=1e-6)


def test_generalized_solver_refuses_a_one_site_chain_rather_than_guessing():
    """Same two-site limit, different consequence. `gs_energy_generalized`
    has no ED fallback, and its outer self-consistent iteration still
    returns a lambda -- the Rayleigh quotient of a state no sweep ever
    touched -- so a one-site chain came back with a silently wrong number
    (-0.3049 against an exact -0.5) rather than an obvious failure. It
    raises now, on the non-Hermitian path too: NH-DMRG escapes ITensor
    v3's short-chain abort because it never calls `dmrg()`, but its own
    sweep is two-site as well, so it does not escape this."""
    for h in ("hermitian", "non_hermitian"):
        sc = spinchain.Spin_Chain(["S=1/2"], itensor_version="python")
        sc.set_hamiltonian(0.8 * sc.Sz[0]
                           + (0.6j if h == "non_hermitian" else 0.6) * sc.Sx[0])
        with pytest.raises(RuntimeError, match="two-site DMRG"):
            sc.gs_energy_generalized(sc.Id)


def test_generalized_solver_still_works_at_two_sites():
    sc = spinchain.Spin_Chain(["S=1/2"] * 2, itensor_version="python")
    sc.set_hamiltonian(sc.Sx[0] * sc.Sx[1] + sc.Sy[0] * sc.Sy[1]
                       + sc.Sz[0] * sc.Sz[1])
    assert sc.gs_energy_generalized(sc.Id) == pytest.approx(
        sc.gs_energy(mode="ED"), abs=1e-6)
