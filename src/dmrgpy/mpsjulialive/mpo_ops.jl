# Algebra of already-built MPOs for mpo.py's MPO class, the julia_live
# counterpart of pyitensor/chain.py's multiply_operators, sum_operators,
# scale_operator, hermitian_operator, trace_operator and
# overlap_aMb_operator.

using LinearAlgebra: tr

function mult_mpo(A,B,maxdim,cutoff)
	# operator product A*B, meaning that B acts first
	return apply(A,B;maxdim=maxdim,cutoff=cutoff)
end

function sum_mpo(A,B,maxdim,cutoff)
	return add(A,B;maxdim=maxdim,cutoff=cutoff)
end

scale_mpo(A,z) = z*A

function dagger_mpo(A)
	# conjugate every tensor and exchange the bra and ket site indices
	return swapprime(dag(A),0=>1;tags="Site")
end

trace_mpo(A) = tr(A) # full trace over the whole Hilbert space

function aMb_mpo(a,A,b)
	# <a|A|b>, exact contraction with no truncation; the stale Link primes
	# an applied operator can leave (see tdvp.jl's tdvp_step) are dropped
	a = noprime(copy(a),"Link")
	b = noprime(copy(b),"Link")
	return inner(prime(a),A,b)
end
