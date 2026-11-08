# Native one-body correlation matrix <Cdag_i C_j> of a spinless fermion
# state, the julia_live counterpart of the session correlation_matrix of
# the other backends (entropytk/correlationentropy.py's
# cpp_correlation_matrix). ITensorMPS's correlation_matrix threads the
# Jordan-Wigner string on Fermion sites itself; it is called on a
# normalised copy, so the result is that of the normalised state.

function correlation_matrix_cdc(psi)
	for s in siteinds(psi)
		if !hastags(s,"Fermion")
			error("correlation_matrix_cdc: every site must be a spinless "*
			      "Fermion site, got an index with tags $(tags(s))")
		end
	end
	phi = normalize(copy(psi))
	return ComplexF64.(correlation_matrix(phi,"Cdag","C"))
end
