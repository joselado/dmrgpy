
# Von Neumann entanglement entropy at a bond, via SVD of the two-site
# tensor there -- the standard MPS way to compute it, and much cheaper
# than reconstructing an explicit density matrix. Mirrors
# pyitensor/chain.py's bond_entropy(wf,b) (itself a port of
# mpscpp3/chain_session.h's Chain::bond_entropy).

function bond_entropy(wf,b)
	psi = copy(wf) # orthogonalize! is in-place; don't mutate the caller's MPS
	orthogonalize!(psi,b)
	twosite = psi[b]*psi[b+1]
	s_b = siteind(psi,b)
	if b>1
		left_link = commonind(psi[b],psi[b-1])
		_,_,_,spec = svd(twosite,left_link,s_b;cutoff=0.0)
	else
		_,_,_,spec = svd(twosite,s_b;cutoff=0.0)
	end
	SvN = 0.0
	for p in eigs(spec)
		if p>1e-12
			SvN += -p*log(p)
		end
	end
	return SvN
end
