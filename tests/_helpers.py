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

    `versions`: which DMRG backend(s) to include (default both). Pass
    `versions=(2,)` for a 2-site chain -- itensor_version=3 crashes hard
    ("LocalOp is default constructed", an ITensor v3 internal check in
    itensor/mps/localop.h) for *any* exactly-2-physical-site chain,
    independent of physics type (confirmed for spin and spinless-fermion
    chains); 3+ sites is unaffected. This is a genuine bug in mpscpp3,
    not a mode.py fallback or a test issue -- see
    test_spin_chain.py/test_fermion_chain.py's dimer tests.
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
