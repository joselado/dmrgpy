import functools


def normalize_terms(terms, ndigits=10):
    """Turn MultiOperator.to_terms() output into an order-independent,
    hashable/sortable form so two term lists can be compared for
    equality regardless of the order terms happen to be generated in.
    """
    return tuple(sorted(
        (round(c.real, ndigits), round(c.imag, ndigits), tuple(tuple(o) for o in ops))
        for c, ops in terms
    ))


def energy_ed_v2_v3(chain, hamiltonian, versions=(2, 3), **kwargs):
    """Ground-state energy of `chain` with Hamiltonian `hamiltonian`,
    computed via ED and via DMRG on both the ITensor v2 and v3 C++
    backends. When a backend isn't compiled, mode.get_mode() silently
    routes that call to ED instead (see mode.py) -- so this still runs
    everywhere, and becomes a real cross-backend regression check
    wherever both extensions are built.

    `chain.set_hamiltonian()` is called fresh before each calculation
    (it also calls `restart()`) because gs_energy()/vev() cache their
    result on the chain (`computed_gs`), which would otherwise just
    return the previous backend's answer instead of recomputing.

    `versions`: which DMRG backend(s) to include (default both). ITensor
    v3's two-site dmrg() used to crash hard ("LocalOp is default
    constructed", an ITensor v3 internal check in itensor/mps/localop.h)
    for any exactly-2-physical-site chain, independent of physics type;
    mode.py's get_mode() now falls back to ED automatically for
    itensor_version==3 with ns<3 (same mechanism as its "extension not
    compiled" fallback), so that no longer crashes -- a v3 request on a
    2-site chain just transparently returns the ED answer instead. The
    `versions` kwarg is kept for tests that want to explicitly restrict
    which backend(s) they exercise (e.g. to keep a "real DMRG on 3+
    sites" test from silently degrading into an ED-only check if the
    chain size ever shrinks), not because including v3 is unsafe.
    """
    chain.set_hamiltonian(hamiltonian)
    e_ed = chain.gs_energy(mode="ED")

    results = [e_ed]
    for version in versions:
        chain.set_hamiltonian(hamiltonian)
        chain.setup_cpp(version)
        results.append(chain.gs_energy(mode="DMRG", **kwargs))

    return tuple(results)


def vev_ed_v2_v3(chain, hamiltonian, observable, versions=(2, 3), **kwargs):
    """Like energy_ed_v2_v3, but for the ground-state expectation value
    of `observable` instead of the ground-state energy itself."""
    chain.set_hamiltonian(hamiltonian)
    v_ed = chain.vev(observable, mode="ED")

    results = [v_ed]
    for version in versions:
        chain.set_hamiltonian(hamiltonian)
        chain.setup_cpp(version)
        results.append(chain.vev(observable, mode="DMRG", **kwargs))

    return tuple(results)


def setup_backend(chain, itensor_version):
    """Switch `chain` onto the requested DMRG backend: itensor_version 2
    or 3 selects the corresponding compiled C++ extension, "python" the
    pure-Python pyitensor engine. Shared by the multi-backend test files
    (each used to hand-roll this same dispatch locally). "julia_live"
    selects the live Julia/ITensors.jl session (mpsjulialive/) -- only
    usable where a Julia toolchain is actually installed, so tests that
    include it should parametrize through julia_live_param() below rather
    than calling this unconditionally."""
    if itensor_version == "python":
        chain.setup_python()
    elif itensor_version == "julia_live":
        chain.setup_julia()
    else:
        chain.setup_cpp(version=itensor_version)


@functools.lru_cache(maxsize=None)
def julia_available():
    """Return (available, reason) for the live Julia backend: available
    is True whenever importing mpsjulialive's juliacall bridge actually
    works, and reason carries the import error otherwise. That import is
    also what boots the Julia session and precompiles ITensors.jl, hence
    the cache -- it is slow the first time and free afterward.

    Mirrors cppext.available(2)/(3) for the compiled C++ extensions, and
    the module-level try/except tests/test_julia_live.py already uses;
    factored out here so the several multi-backend test files that
    parametrize over julia_live share one probe instead of each paying
    for (and hand-rolling) their own."""
    try:
        from dmrgpy.mpsjulialive import juliasession  # noqa: F401
        return True, ""
    except Exception as e:  # pragma: no cover - environment dependent
        return False, str(e)


def julia_live_param():
    """pytest.param("julia_live", ...) that skips itself when no Julia
    toolchain is available -- the julia_live counterpart of the
    pytest.param(3, marks=skipif(not cppext.available(3))) idiom used
    throughout this suite."""
    import pytest
    ok, reason = julia_available()
    return pytest.param("julia_live", id="julia_live",
            marks=pytest.mark.skipif(not ok,
                reason="requires a working juliacall/Julia toolchain: %s"
                        % reason))
