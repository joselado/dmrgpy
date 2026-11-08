# Non-Hermitian KPM (port of NHKPM.jl): the biorthogonal Chebyshev moments
# of nonhermitian/kpm.py, computed from a coupled forward/adjoint
# recursion built from hs=(z*Id-H)/E_max and its dagger hsdag. Line by
# line the same recursion as pyitensor/chain.py's Chain.nhkpm_moments and
# mpscpp3/chain_session.h's Chain::nhkpm_moments (see the latter for the
# full algorithm comment). Uses mpsalgebra.jl's apply_op and summps, so it
# must load after it. Returns mu, length n, mu[k] = <wfb|vn_{2k-1}> in the
# 1-based notation of the reference.

function nhkpm_moments(hs,hsdag,wfa,wfb,n,maxdim,cutoff)
	if n < 1
		error("nhkpm_moments: n must be >= 1")
	end
	v = 1.0*wfa
	alpha_prev2 = apply_op(hsdag,v,maxdim,cutoff) # alpha[1]
	alpha_prev1 = summps(2.0*apply_op(hs,alpha_prev2,maxdim,cutoff),
			(-1.0)*v,maxdim,cutoff) # alpha[2]
	vn_prev2 = 1.0*v # vn[1]
	vn_prev1 = 2.0*apply_op(hs,vn_prev2,maxdim,cutoff) # vn[2]
	mu = ComplexF64[inner(wfb,vn_prev2)]
	for k=2:n
		alpha_x = summps(2.0*apply_op(hsdag,alpha_prev1,maxdim,cutoff),
				(-1.0)*alpha_prev2,maxdim,cutoff)
		alpha_y = summps(2.0*apply_op(hs,alpha_x,maxdim,cutoff),
				(-1.0)*alpha_prev1,maxdim,cutoff)
		t = summps(2.0*apply_op(hsdag,vn_prev1,maxdim,cutoff),
				(-1.0)*vn_prev2,maxdim,cutoff)
		vn_x = summps(2.0*alpha_prev1,t,maxdim,cutoff)
		vn_y = summps(2.0*apply_op(hs,vn_x,maxdim,cutoff),
				(-1.0)*vn_prev1,maxdim,cutoff)
		push!(mu,inner(wfb,vn_x))
		alpha_prev2,alpha_prev1 = alpha_x,alpha_y
		vn_prev2,vn_prev1 = vn_x,vn_y
	end
	return mu
end
