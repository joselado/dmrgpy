
# Single-site reduced density matrix, via explicit progressive tracing of
# every site to the right of `site` -- mirrors pyitensor/chain.py's own
# reduced_dm(wf,site) (itself a port of mpscpp3/chain_session.h's
# Chain::reduced_dm).
#
# The state is divided by its norm, sqrt(<psi|psi>). It used to be divided
# by <psi|psi> itself, a quirk ported verbatim from the other backends and
# commented as a no-op because a solved state is unit norm; set_gs(c*s)
# and gs_energy(wf0=c*s, reconverge=False) hand over a caller's state
# unswept, and the matrix came back as rho/c^2, trace 0.25 at c=2
# (2026-09-25 hole hunt 25b, finding 11).
#
# At the last site there is no right link to prime: psi[site+1] is out of
# range there (a BoundsError, finding 12 of the same record, the Julia half
# of 2026-08 finding 16), and priming the physical index alone is the same
# contraction, since orthogonalize!(psi,site) has already put everything to
# the left into psi[site] -- the guard pyitensor/chain.py got for that
# finding, line for line.

function reduced_dm(wf,site)
	psi = copy(wf)
	nrm2 = real(inner(psi,psi))
	psi = psi*(1.0/sqrt(nrm2))
	orthogonalize!(psi,site)
	s = siteind(psi,site)
	if site < length(psi)
		ir = commonind(psi[site],psi[site+1])
		rho = psi[site]*dag(prime(psi[site],s,ir))
	else
		rho = psi[site]*dag(prime(psi[site],s))
	end
	for k=site+1:length(psi)
		rho = rho*psi[k]
		rho = rho*dag(prime(psi[k],"Link"))
	end
	return Array(rho,s,prime(s))
end
