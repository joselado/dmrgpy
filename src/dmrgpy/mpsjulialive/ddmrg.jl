# Dynamical DMRG (Jeckelmann, Phys. Rev. B 66, 045114): the correction
# vector of cvm_solver="variational" on itensor_version="julia_live".
#
# A line-by-line port of mpscpp3's Chain::ddmrg_correction_vector (and of
# pyitensor/ddmrg.py): a two-site sweep minimizing
#
#     W(x) = <x|M|x> - 2 Re<b|x>,   M = K^2 + eta^2,   K = H - E0 - omega,
#
# with b = -eta B|GS>. Each two-site block is solved by conjugate gradient
# on the projected system (Hermitian positive definite, nothing truncated
# inside it), then split by SVD, so W never increases from sweep to sweep.
#
# Two differences from the v3 code, both exact:
#   * M is never built as an MPO. In the mixed canonical gauge of the sweep
#     the projection of eta^2*Id onto the two-site block is eta^2 times the
#     identity, so the local matvec adds eta^2*theta, and the global W adds
#     eta^2*<x|x>. Only K^2 is an MPO, so no MPO sum is needed.
#   * The caller (mpsjulialive/ddmrg.py) hands K, eta and b already in the
#     units of the solver scale rule, so this file has no scale of its own.
#
# Uses mpsalgebra.jl's apply_op, so it is loaded after it.

using ITensors
using ITensorMPS
using LinearAlgebra


# Conjugate gradient on one two-site block, M_eff theta = rhs. Same
# tolerance (1e-10 of ||rhs||) and iteration cap (200) as v3 and python.
function ddmrg_local_cg(mv, rhs, x1; tol = 1e-10, maxiter = 200)
	bn = norm(rhs)
	bn == 0.0 && return 0.0 * rhs
	r = rhs - mv(x1)
	p = r
	rs = norm(r)^2
	for it in 1:maxiter
		sqrt(rs) <= tol * bn && break
		Ap = mv(p)
		den = real(scalar(dag(p) * Ap))
		den == 0.0 && break
		alpha = rs / den
		x1 = x1 + alpha * p
		r = r - alpha * Ap
		rsn = norm(r)^2
		p = r + (rsn / rs) * p
		rs = rsn
	end
	return x1
end


# K^2 as an MPO: the exact product, then compressed. The compression only
# removes redundant channels, hence the tight build cutoff.
function ddmrg_square_mpo(K, mpomaxm)
	K2 = apply(K, K; alg = "naive", truncate = false)
	truncate!(K2; cutoff = 1e-14, maxdim = mpomaxm)
	return K2
end


# The correction vector sweep. Returns (W, x), with W the functional
# W = <x|K^2|x> + eta^2 <x|x> - 2 Re<b|x> of the returned state x.
#   K       MPO of (H - E0 - omega), in solver units
#   S2      MPO of B
#   psi0    the ground state
#   eta     the broadening, in solver units
#   bscale  b = bscale * B|GS>
#   x0      a warm start (the previous frequency's x), or nothing
function ddmrg_correction_vector(K, S2, psi0, eta, bscale, maxdim, nsweeps,
		cutoff, mpomaxm, x0 = nothing)
	eta2 = eta * eta
	psi = noprime(copy(psi0), "Link") # a state built by operator algebra
	                                  # can carry a stale Link prime
	K2 = ddmrg_square_mpo(K, mpomaxm)
	b = apply_op(S2, psi, maxdim, cutoff)
	b = noprime(b, "Link")
	b = bscale * b
	# x shares no Link index with b: start from a fresh copy
	x = x0 === nothing ? copy(b) : noprime(copy(x0), "Link")
	x = sim(linkinds, x)
	n = length(x)
	orthogonalize!(x, 1)
	one_ = ITensor(1.0)
	# K^2 environments (bra x primed, ket x) and overlap environments <x|b>
	LM = Dict{Int, ITensor}(0 => one_)
	RM = Dict{Int, ITensor}(n + 1 => one_)
	Lo = Dict{Int, ITensor}(0 => one_)
	Ro = Dict{Int, ITensor}(n + 1 => one_)
	envM(prev, i) = prev * x[i] * K2[i] * dag(prime(x[i]))
	envO(prev, i) = prev * b[i] * dag(prime(x[i], "Link"))
	function local_update(i, right)
		theta = x[i] * x[i + 1]
		L = LM[i - 1]
		R = RM[i + 2]
		mv(v) = noprime(L * v * K2[i] * K2[i + 1] * R) + eta2 * v
		rhs = Lo[i - 1] * b[i] * b[i + 1] * Ro[i + 2]
		rhs = noprime(rhs)
		theta = ddmrg_local_cg(mv, rhs, theta)
		if right
			li = commoninds(theta, x[i])
			U, S, V = svd(theta, li; cutoff = cutoff, maxdim = maxdim,
				lefttags = "Link,l=$i", righttags = "Link,l=$i")
			x[i] = U
			x[i + 1] = S * V
		else
			li = commoninds(theta, x[i + 1])
			U, S, V = svd(theta, li; cutoff = cutoff, maxdim = maxdim,
				lefttags = "Link,l=$i", righttags = "Link,l=$i")
			x[i + 1] = U
			x[i] = V * S
		end
	end
	for sw in 1:nsweeps
		for i in n:-1:3
			RM[i] = envM(RM[i + 1], i)
			Ro[i] = envO(Ro[i + 1], i)
		end
		for i in 1:(n - 1)
			local_update(i, true)
			if i < n - 1
				LM[i] = envM(LM[i - 1], i)
				Lo[i] = envO(Lo[i - 1], i)
			end
		end
		for i in (n - 1):-1:1
			local_update(i, false)
			if i > 1
				RM[i + 1] = envM(RM[i + 2], i + 1)
				Ro[i + 1] = envO(Ro[i + 2], i + 1)
			end
		end
	end
	W = real(inner(x', K2, x)) + eta2 * real(inner(x, x)) -
		2.0 * real(inner(b, x))
	return W, x
end
