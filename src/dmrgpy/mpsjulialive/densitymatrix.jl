
# Single-site reduced density matrix, via explicit progressive tracing of
# every site to the right of `site` -- mirrors pyitensor/chain.py's own
# reduced_dm(wf,site) (itself a port of mpscpp3/chain_session.h's
# Chain::reduced_dm), including its "divide by the norm squared, not its
# square root" quirk (preserved for cross-backend consistency -- in
# practice a no-op, since wf is essentially always already unit-norm
# coming out of dmrg()).

function reduced_dm(wf,site)
	psi = copy(wf)
	nrm2 = real(inner(psi,psi))
	psi = psi*(1.0/nrm2)
	orthogonalize!(psi,site)
	ir = commonind(psi[site],psi[site+1])
	s = siteind(psi,site)
	rho = psi[site]*dag(prime(psi[site],s,ir))
	for k=site+1:length(psi)
		rho = rho*psi[k]
		rho = rho*dag(prime(psi[k],"Link"))
	end
	return Array(rho,s,prime(s))
end
