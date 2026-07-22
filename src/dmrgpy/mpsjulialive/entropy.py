def bond_entropy(psi,b):
    """Von Neumann entanglement entropy at bond b (1-indexed, between
    sites b and b+1), on the Julia backend -- mpsjulialive/entropy.jl's
    bond_entropy, a native SVD-based computation."""
    from .juliasession import Main as Mainjl
    return Mainjl.bond_entropy(psi.jlmps,b)
