def normalize_terms(terms, ndigits=10):
    """Turn MultiOperator.to_terms() output into an order-independent,
    hashable/sortable form so two term lists can be compared for
    equality regardless of the order terms happen to be generated in.
    """
    return tuple(sorted(
        (round(c.real, ndigits), round(c.imag, ndigits), tuple(tuple(o) for o in ops))
        for c, ops in terms
    ))


def energy_ed_v2_v3(chain, hamiltonian, **kwargs):
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
    """
    chain.set_hamiltonian(hamiltonian)
    e_ed = chain.gs_energy(mode="ED")

    chain.set_hamiltonian(hamiltonian)
    chain.setup_cpp(2)
    e_v2 = chain.gs_energy(mode="DMRG", **kwargs)

    chain.set_hamiltonian(hamiltonian)
    chain.setup_cpp(3)
    e_v3 = chain.gs_energy(mode="DMRG", **kwargs)

    return e_ed, e_v2, e_v3


def vev_ed_v2_v3(chain, hamiltonian, observable, **kwargs):
    """Like energy_ed_v2_v3, but for the ground-state expectation value
    of `observable` instead of the ground-state energy itself."""
    chain.set_hamiltonian(hamiltonian)
    v_ed = chain.vev(observable, mode="ED")

    chain.set_hamiltonian(hamiltonian)
    chain.setup_cpp(2)
    v_v2 = chain.vev(observable, mode="DMRG", **kwargs)

    chain.set_hamiltonian(hamiltonian)
    chain.setup_cpp(3)
    v_v3 = chain.vev(observable, mode="DMRG", **kwargs)

    return v_ed, v_v2, v_v3
