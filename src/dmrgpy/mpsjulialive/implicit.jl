# Implicit imaginary-time steps and linear solves of an MPS by two-site
# Galerkin sweeps, on itensor_version="julia_live". A port of
# pyitensor/implicit.py (arXiv:2606.02930), whose module docstring has the
# math; the functions here mirror its functions one by one.
#
# One sweep (implicit_sweep!) moves x towards the solution of
#
#     (a0 + c0 H) x = b0 psi + d0 H psi,
#
# projected at each bond (i,i+1) on the tangent space of x. With x in
# mixed canonical gauge the projected identity is the identity, so the
# local problem is (a0 + c0 H_eff) theta = b0 P[psi] + d0 P[H psi], with
# P[psi] from the overlap environments <x|psi> and P[H psi] from the mixed
# environments <x|H|psi>: no MPO of a0 + c0 H and no MPS of the right-hand
# side is built. The local problem is solved by restarted GMRES
# (KrylovKit) and theta is split by SVD.
#
# The environments are written by hand with ITensors primitives, as in
# ddmrg.jl, and not through ProjMPO or ITensorMPS's alternating_update,
# whose internals change between versions. x carries link indices of its
# own (sim), so an environment between x and psi contracts only the site
# indices. Dense sites only: the Python side refuses a conserved sector.

using ITensors
using ITensorMPS
using LinearAlgebra

# KrylovKit is declared in juliapkg.json and loaded by its UUID, as in
# vumps.jl, which returns the copy ITensorMPS has loaded already; the
# fallback is the binding inside ITensorMPS that cbe.jl uses.
const _IMPLICIT_KK = try
	Base.require(Base.PkgId(
		Base.UUID("0b1a1467-8014-51b9-945f-bf0ae24f4b77"), "KrylovKit"))
catch
	ITensorMPS.KrylovKit
end


# A state without the stale prime level that operator algebra leaves on
# its Link indices (tdvp.jl's note), copied, so that the caller's MPS is
# never touched.
implicit_clean(psi) = noprime(copy(psi), "Link")


# (a0 + c0 heff) theta = rhs by restarted GMRES from theta0, stopped at
# ||rhs - A theta|| <= tol ||rhs||. As in pyitensor's gmres, krylov is the
# Krylov dimension of one restart and maxiter counts matrix-vector
# products; KrylovKit's maxiter counts restarts, each of one residual and
# up to krylov products.
function implicit_local_solve(heff, rhs, theta0, a0, c0, tol, krylov, maxiter)
	bn = norm(rhs)
	bn == 0.0 && return 0.0 * rhs
	restarts = max(1, cld(maxiter - 1, krylov + 1))
	theta, _ = _IMPLICIT_KK.linsolve(heff, rhs, theta0, a0, c0;
		ishermitian = false, krylovdim = krylov, maxiter = restarts,
		rtol = tol, atol = 0.0, verbosity = 0)
	return theta
end


# One two-site sweep of x, a right and a left half sweep, in place (see
# the header). x has its orthogonality centre on site 1 on entry and on
# exit. The mixed environments are built only when d0 is not zero.
function implicit_sweep!(x, psi, H, a0, c0, b0, d0, maxdim, cutoff, tol,
		krylov, maxiter)
	n = length(x)
	mixed = d0 != 0
	one_ = ITensor(1.0)
	envH(prev, i) = prev * x[i] * H[i] * dag(prime(x[i]))   # <x|H|x>
	envM(prev, i) = prev * psi[i] * H[i] * dag(prime(x[i])) # <x|H|psi>
	envO(prev, i) = prev * psi[i] * dag(x[i])               # <x|psi>
	LH = Dict{Int, ITensor}(0 => one_)
	LM = Dict{Int, ITensor}(0 => one_)
	LO = Dict{Int, ITensor}(0 => one_)
	RH = Dict{Int, ITensor}(n + 1 => one_)
	RM = Dict{Int, ITensor}(n + 1 => one_)
	RO = Dict{Int, ITensor}(n + 1 => one_)
	function extend!(H_, M_, O_, prevkey, i)
		H_[i] = envH(H_[prevkey], i)
		mixed && (M_[i] = envM(M_[prevkey], i))
		O_[i] = envO(O_[prevkey], i)
	end
	# ortho = "left" keeps site i left-orthogonal and moves the centre to
	# i+1 (pyitensor's direction "right"), "right" the mirror
	function update!(i, ortho)
		L = LH[i - 1]
		R = RH[i + 2]
		heff(v) = noprime(L * v * H[i] * H[i + 1] * R)
		rhs = b0 * (LO[i - 1] * psi[i] * psi[i + 1] * RO[i + 2])
		if mixed
			rhs = rhs + d0 * noprime(LM[i - 1] * psi[i] * H[i] * psi[i + 1] *
				H[i + 1] * RM[i + 2])
		end
		theta = implicit_local_solve(heff, rhs, x[i] * x[i + 1], a0, c0, tol,
			krylov, maxiter)
		replacebond!(x, i, theta; maxdim = maxdim, cutoff = cutoff,
			ortho = ortho, which_decomp = "svd")
	end
	for i in n:-1:3
		extend!(RH, RM, RO, i + 1, i)
	end
	for i in 1:(n - 1)
		update!(i, "left")
		i < n - 1 && extend!(LH, LM, LO, i - 1, i)
	end
	for i in (n - 1):-1:1
		update!(i, "right")
		i > 1 && extend!(RH, RM, RO, i + 2, i + 1)
	end
	return x
end


# x = (1 + beta(H - shift))^-1 (1 - beta(H - shift)) psi by nsweeps
# sweeps from psi itself (fresh links). Not normalized, orthogonality
# centre on site 1.
function implicit_factor(psi, H, beta, shift, maxdim, cutoff, nsweeps, tol,
		krylov, maxiter)
	a0 = 1.0 - beta * shift # theta coefficient of the local operator
	b0 = 1.0 + beta * shift # psi coefficient of the right-hand side
	x = sim(linkinds, psi)
	orthogonalize!(x, 1)
	for _ in 1:nsweeps
		implicit_sweep!(x, psi, H, a0, beta, b0, -beta, maxdim, cutoff, tol,
			krylov, maxiter)
	end
	return x
end


# One imaginary-time step: the factors of betas in turn, then the state
# normalized. shift_energy takes the shift from <psi|H|psi>/<psi|psi> of
# the input state. Mirrors pyitensor's implicit_step and mpscpp3's
# Chain::implicit_step.
function implicit_step(psi, H, betas, shift, shift_energy, maxdim, cutoff,
		nsweeps, tol, krylov, maxiter)
	(nsweeps >= 1 && krylov >= 1 && maxiter >= 1) ||
		error("implicit_step: nsweeps, krylov and maxiter must be >= 1")
	psi = implicit_clean(psi)
	length(psi) >= 2 || error("implicit_step needs at least 2 sites")
	s = shift_energy ? real(inner(psi', H, psi)) / real(inner(psi, psi)) :
		Float64(shift)
	x = psi
	for b in betas
		x = implicit_factor(x, H, ComplexF64(b), s, maxdim, cutoff, nsweeps,
			tol, krylov, maxiter)
	end
	nrm = norm(x)
	nrm > 0.0 || error("implicit_step: the state vanished")
	return x / nrm
end


# ||A x - psi||, with A x and the difference formed exactly. The norm is
# read off the orthogonality centre of the difference, since <r|r> from
# the full contraction is a cancellation that floors near 1e-8 ||psi||
# (pyitensor's _residual).
function implicit_residual(A, x, psi)
	Ax = apply(A, x; alg = "naive", truncate = false)
	r = add(Ax, (-1.0) * psi; alg = "directsum")
	orthogonalize!(r, 1)
	return norm(r[1])
end


# x solving A x = psi for an MPO A by two-site sweeps, A_eff theta = P[psi]
# at each bond. The sweeps stop once ||A x - psi|| <= tol (absolute), after
# nsweeps sweeps, or once a sweep no longer lowers the residual by a tenth.
# Returns the iterate with the smallest residual, not normalized, and the
# residual history. Mirrors pyitensor's linear_solve, with its defaults.
function implicit_linear_solve(A, psi, tol, nsweeps, maxdim, cutoff;
		local_tol = 1e-10, krylov = 30, maxiter = 300)
	psi = implicit_clean(psi)
	length(psi) >= 2 || error("implicit_linear_solve needs at least 2 sites")
	x = sim(linkinds, psi)
	orthogonalize!(x, 1)
	history = Float64[implicit_residual(A, x, psi)]
	best, best_res = copy(x), history[1]
	for _ in 1:nsweeps
		best_res <= tol && break
		implicit_sweep!(x, psi, A, 0.0, 1.0, 1.0, 0.0, maxdim, cutoff,
			local_tol, krylov, maxiter)
		res = implicit_residual(A, x, psi)
		push!(history, res)
		improved = res < 0.9 * best_res
		if res < best_res
			best, best_res = copy(x), res
		end
		improved || break
	end
	return best, history
end
