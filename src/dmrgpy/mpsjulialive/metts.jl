
# METTS (Minimally Entangled Typical Thermal States; E.M. Stoudenmire and
# S.R. White, "Minimally entangled typical thermal state algorithms", New
# J. Phys. 12, 055026 (2010), arXiv:1002.1305), run entirely in Julia (one
# call per full nwarmup+nsamples Markov chain, not one Python<->Julia round
# trip per sample -- same design as kpm.jl's moment recursion and tdvp.jl's
# quench_tdvp/evolve_and_measure_tdvp). A value-level port of
# pyitensor/metts.py's algorithm (see that module's own docstring for the
# full derivation of why this Markov chain's stationary distribution
# already carries the correct Boltzmann weight, so a plain sample average
# converges to the thermal average with no importance reweighting) --
# mpscpp3/chain_session.h's Chain::metts_vev is the same port onto real
# ITensor v3 via its own diagHermitian()/setElt() primitives; this file
# reaches the same result by extracting plain Julia arrays at each site
# (LinearAlgebra.eigen on a Hermitian matrix, explicit array bookkeeping)
# rather than through ITensor-level tensor primitives, since that is the
# safest direct translation of the already-validated pyitensor reference
# rather than a fresh ITensor-native derivation.
#
# Reuses tdvp.jl's tdvp_step() unchanged for the imaginary-time evolution
# step: that function's own -im*dt convention already generalizes to a
# purely real dt (its usual real-time case) or -- what is wanted here -- a
# purely imaginary dt, since -im*(-im*dtau) = -dtau gives genuine decay
# rather than rotation (same reasoning as tdz.jl's
# advance_complex_time_step, ../timedependent.py:89). No new evolution
# primitive is needed; METTS only adds the sampling/collapse machinery
# below.

using ITensors
using ITensorMPS
using Random
using LinearAlgebra


# The (dim,dim) matrix of single-site operator `opname` at site `s`, in
# the convention that plain matrix-vector multiplication reproduces
# op(opname,s) acting on a bare physical-index vector -- obtained by
# reading off the ITensor's own (primed-out, unprimed-in) index
# convention via Array(), rather than re-deriving sites.jl's operator
# tables algebraically here.
function metts_site_operator_matrix(s, opname)
	O = op(opname, s)
	return Array(O, s', s)
end


# (evals, U) for the named single-site Hermitian operator at site s: evals
# ascending, U's columns the corresponding eigenvectors expressed in the
# site's native (computational) basis -- every physical operator dmrgpy's
# site tables define (Sz, Sx, N, ...) is Hermitian, so eigen(Hermitian(.))
# is exact (not merely a numerically-stable approximation).
function metts_eigenbasis(s, opname)
	M = metts_site_operator_matrix(s, opname)
	F = eigen(Hermitian(M))
	return F.values, F.vectors
end


# {(opname,i) => (evals,U)} for every (opname,site) pair a run will ever
# need, computed once up front -- mirrors pyitensor/metts.py's
# _build_eigcache()/mpscpp3's metts_build_eigcache(): a given (opname,i)'s
# eigenbasis never changes across the whole sampling run, so this avoids
# tens of thousands of redundant eigendecompositions across
# nwarmup+nsamples iterations.
function metts_build_eigcache(sites, basis_ops)
	cache = Dict{Tuple{String,Int},Tuple{Vector{Float64},Matrix{ComplexF64}}}()
	n = length(sites)
	for opname in basis_ops
		for i in 1:n
			key = (opname, i)
			if haskey(cache, key)
				continue
			end
			cache[key] = metts_eigenbasis(sites[i], opname)
		end
	end
	return cache
end


# A classical product state (CPS: bond dimension 1 throughout) from a
# per-site vector of amplitudes. Fresh dim-1 Link indices are built
# between neighboring sites so this looks, to tdvp_step/inner/orthogonalize!
# and every other downstream routine, like any other MPS -- mirrors
# pyitensor/metts.py's product_state() and mpscpp3's build_product_state().
function metts_product_state(sites, vectors)
	n = length(sites)
	if n == 1
		return MPS([ITensor(vectors[1], sites[1])])
	end
	links = [Index(1, "Link,l=$k") for k in 1:(n-1)]
	tensors = ITensor[]
	for i in 1:n
		s = sites[i]
		v = vectors[i]
		if i == 1
			T = ITensor(reshape(v, :, 1), s, links[1])
		elseif i == n
			T = ITensor(reshape(v, 1, :), links[n-1], s)
		else
			T = ITensor(reshape(v, 1, :, 1), links[i-1], s, links[i])
		end
		push!(tensors, T)
	end
	return MPS(tensors)
end


# A random CPS: at each site, uniformly picks one eigenstate of the named
# single-site Hermitian operator. Valid seed for the METTS Markov chain
# regardless of choice -- the chain's stationary distribution doesn't
# depend on the starting CPS (paper, Sec. 2), only how many warmup steps
# are needed to reach it.
function metts_random_cps(sites, opname, eigcache, rng)
	n = length(sites)
	vectors = Vector{Vector{ComplexF64}}(undef, n)
	for i in 1:n
		_, U = eigcache[(opname, i)]
		d = size(U, 2)
		k = rand(rng, 1:d)
		vectors[i] = U[:, k]
	end
	return metts_product_state(sites, vectors)
end


# Inverse-CDF sample from a discrete distribution p (already normalized to
# sum to 1) -- the Julia-stdlib-only equivalent of numpy's
# rng.choice(d,p=p) / std::discrete_distribution, avoiding a Distributions.jl
# dependency (not declared in ../juliapkg.json) for this one use. Falls
# back to the last outcome on roundoff (u landing fractionally above the
# cumulative sum's final value of ~1), same convention as re-normalizing p
# in the callers below.
function metts_sample_categorical(rng, p)
	u = rand(rng)
	c = 0.0
	for k in 1:length(p)
		c += p[k]
		if u <= c
			return k
		end
	end
	return length(p)
end


# Sequential ("perfect sampling") collapse of MPS wf onto a new CPS built
# from eigenstates of the single-site Hermitian operator opname -- direct
# port of pyitensor/metts.py's collapse_to_cps() (see its own docstring
# for the full derivation), the paper's "quantum measurement" collapse
# step (Sec. 2.2) generalized to the whole chain via one left-to-right
# sweep.
#
# orthogonalize(wf,1) first makes every site i>1 right-orthonormal;
# contracting the running "collapsed-so-far" amplitude L into site i's
# tensor and rotating into the eigenbasis (U' * mat) then gives exactly
# the marginal probability of outcome k at site i conditioned on the
# outcomes already sampled to its left -- no explicit partial trace
# needed, right-orthonormality already performs it.
function metts_collapse_to_cps(wf, sites, opname, eigcache, rng)
	n = length(sites)
	psi = orthogonalize(wf, 1)
	vectors = Vector{Vector{ComplexF64}}(undef, n)
	L = nothing  # running collapsed-prefix amplitude (ITensor over the current left link), or nothing at site 1
	for i in 1:n
		evals, U = eigcache[(opname, i)]
		s = sites[i]
		d = length(evals)
		T = psi[i]
		if L !== nothing
			T = L * T
		end
		right_link = i < n ? commonind(psi[i], psi[i+1]) : nothing
		if right_link !== nothing
			mat = Array(T, s, right_link)
		else
			mat = reshape(Array(T, s), d, 1)
		end
		rot = U' * mat  # amplitude of each eigenbasis outcome, per remaining right-link value
		probs_raw = vec(sum(abs2.(rot), dims=2))
		total = sum(probs_raw)
		p = probs_raw ./ total
		p = p ./ sum(p)  # re-normalize away roundoff
		k = metts_sample_categorical(rng, p)
		vectors[i] = U[:, k]
		if right_link !== nothing
			Lvec = rot[k, :] ./ sqrt(probs_raw[k])
			L = ITensor(Lvec, right_link)
		end
	end
	return metts_product_state(sites, vectors)
end


# Evolve psi by e^{-dbeta_half*H}, via nsteps of imaginary-time TDVP
# (tdvp_step's own -im*dt convention with dt purely imaginary; see this
# file's header), renormalizing after every step since imaginary-time
# evolution isn't norm-preserving (unlike tdvp_step's usual real-time use).
function metts_imaginary_time_evolve(psi, H, dbeta_half, nsteps, cutoff, maxdim)
	dt = -im * (dbeta_half / nsteps)
	for _ in 1:nsteps
		psi = tdvp_step(H, psi, dt, cutoff, maxdim)
		psi = psi / sqrt(abs(inner(psi, psi)))
	end
	return psi
end


# Thermal <Op> for every op in opmpos, via METTS sampling -- the
# julia_live counterpart of pyitensor/metts.py's metts_thermal_average()/
# mpscpp3/chain_session.h's Chain::metts_vev(), see either for the
# algorithm's own derivation. Runs the whole nwarmup+nsamples Markov chain
# in this one call (see file header). opmpos: a Julia/Python-list of MPOs
# (any object usable as tdvp_step's H / inner()'s middle operand);
# basis_ops: a list of single-site operator name strings, cycled one per
# sample for ergodicity (paper, Sec. 3.3).
#
# Returns (means,stderrs): means[k] is the sample mean of opmpos[k],
# stderrs[k] its naive iid standard error over the nsamples retained
# samples -- likely optimistic since consecutive METTS samples are
# Markov-correlated (see pyitensor/metts.py's own docstring).
function metts_vev(H, sites, opmpos, T, nsamples, nwarmup, dbeta_half_step,
		basis_ops, seed, cutoff, maxdim)
	beta_half = 1.0 / (2.0 * T)
	nsteps = max(1, Int(ceil(beta_half / dbeta_half_step)))
	rng = MersenneTwister(seed)
	nbasis = length(basis_ops)
	eigcache = metts_build_eigcache(sites, basis_ops)
	cps = metts_random_cps(sites, basis_ops[1], eigcache, rng)

	nops = length(opmpos)
	samples = [ComplexF64[] for _ in 1:nops]
	total_iters = nwarmup + nsamples
	for it in 1:total_iters
		phi = metts_imaginary_time_evolve(cps, H, beta_half, nsteps, cutoff, maxdim)
		if it > nwarmup
			for k in 1:nops
				push!(samples[k], inner(phi, opmpos[k], phi))
			end
		end
		basis = basis_ops[((it - 1) % nbasis) + 1]
		cps = metts_collapse_to_cps(phi, sites, basis, eigcache, rng)
	end

	means = ComplexF64[]
	stderrs = Float64[]
	for k in 1:nops
		vals = samples[k]
		m = sum(vals) / length(vals)
		push!(means, m)
		if length(vals) > 1
			var = sum(abs2.(vals .- m)) / (length(vals) - 1)
			push!(stderrs, sqrt(var / length(vals)))
		else
			push!(stderrs, 0.0)
		end
	end
	return means, stderrs
end
