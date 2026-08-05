
# Generalized-eigenvalue DMRG on the live Julia backend: the smallest
# lambda solving H|psi> = lambda*A|psi> for a Hermitian positive-definite
# metric MPO A, plus its non-Hermitian (complex-lambda, biorthogonal)
# counterpart.
#
# Both are ports of the same self-consistent Lagrange-multiplier
# algorithm already implemented for the pure-Python backend
# (pyitensor/dmrg.py's dmrg_generalized(), pyitensor/nhdmrg.py's
# nhdmrg_generalized()) and for the compiled ITensor v3 one
# (mpscpp3/chain_session.h's Chain::gs_energy_generalized/
# Chain::nhdmrg_generalized) -- see those docstrings for the derivation.
# In short: minimizing <psi|H|psi> under the *metric* normalization
# <psi|A|psi>=1 has stationarity condition (H-lambda*A)|psi> = mu|psi>,
# i.e. the ordinary (plain-normalized) eigenproblem of the shifted
# operator H-lambda*A, which is exactly what one standard two-site DMRG
# sweep already solves. So each outer iteration (i) rebuilds
# Heff = H - lambda*A from the current lambda, (ii) runs one ordinary
# sweep against it, and (iii) resets lambda to the freshly-swept state's
# generalized Rayleigh quotient <psi|H|psi>/<psi|A|psi>. At convergence
# mu -> 0 and the stationarity condition is precisely
# H|psi> = lambda*A|psi>.
#
# Written as one Julia function per method rather than a Python-side loop
# calling per-sweep Julia primitives, for the same reason kpm.jl and
# tdvp.jl are: the whole iteration stays inside one Julia call, with no
# Python<->Julia round trip (nor MPO/MPS marshalling) per outer sweep.
# Sweeps construction and stdout silencing are get_gs.jl's shared
# make_sweeps()/run_quiet() (loaded before this file, see
# juliasession.py's initialize()).

using ITensors
using ITensorMPS
using ITensorNHDMRG


# H and A are ordinary Hamiltonian-like MPOs of modest bond dimension, so
# their shifted sum is formed *exactly* (cutoff=0.0, no maxdim) rather
# than truncated: it costs little, it never accumulates across outer
# iterations (Heff is always rebuilt from the original H and A, never
# from the previous Heff), and truncating the shift term would silently
# change the operator whose spectrum the sweep is targeting. Mirrors
# dmrg_generalized()'s own exact mps_sum(H, A*(-lam)).
function shifted_mpo(H, A, lam)
	return add(H, (-lam) * A; cutoff = 0.0)
end


# Guard shared by both solvers: <psi|A|psi> (or its biorthogonal
# counterpart <psi_L|A|psi_R>) collapsing to ~0 mid-iteration means A
# isn't actually positive definite for this problem, or the sweep drove
# the state into A's near-null-space. Dividing through it would poison
# every later outer iteration with Inf/NaN instead of failing loudly.
# Written as !(abs(a)>tol) rather than abs(a)<tol so a NaN -- which fails
# *both* comparisons -- is caught too, same as the Python/C++ ports'
# identical guards.
function check_metric_norm(a, what)
	if !(abs(a) > 1e-14)
		error(what * ": the metric expectation value collapsed to ~0 " *
		      "(or NaN) mid-iteration (A may not be positive definite, " *
		      "or the sweep drove the state toward A's near-null-space) " *
		      "-- cannot form a meaningful generalized Rayleigh quotient")
	end
end


"""
    get_gs_generalized(H, A, psi0; nsweeps, cutoff, maxm, lam0)

Hermitian generalized eigenproblem. Returns `(lam, psi)`.

`lam0` is the starting lambda estimate; `NaN` means "unset" and defaults
to `psi0`'s own generalized Rayleigh quotient (two cheap `inner`s, giving
the first outer sweep a head start over starting from lambda=0). NaN is
the sentinel rather than `nothing` so this signature stays interchangeable
with `Chain::gs_energy_generalized`'s own pybind11 binding, which has no
`None` equivalent either.
"""
function get_gs_generalized(H, A, psi0; nsweeps = 10, cutoff = 1e-8,
		maxm = 80, lam0 = NaN)
	psi = psi0
	lam = lam0
	if isnan(lam)
		a0 = real(inner(psi', A, psi))
		h0 = real(inner(psi', H, psi))
		lam = abs(a0) > 1e-14 ? h0 / a0 : 0.0
	end
	# One outer iteration == one ordinary DMRG sweep, so the schedule is
	# a fresh single-sweep Sweeps rebuilt against each new Heff rather
	# than one nsweeps-long Sweeps handed to a single dmrg() call.
	sweeps = make_sweeps(1, maxm, cutoff)
	for i = 1:nsweeps
		Heff = shifted_mpo(H, A, lam)
		_e, psi = run_quiet() do
			dmrg(Heff, psi, sweeps)
		end
		a_psi = real(inner(psi', A, psi))
		check_metric_norm(a_psi, "get_gs_generalized")
		lam = real(inner(psi', H, psi)) / a_psi
	end
	return lam, psi
end


# <psil|X|psir> for an ITensorNHDMRG (wfl,wfr) pair, in dmrgpy's own
# left-vector convention (psil = conj(wfl), see nhdmrg.jl's header) --
# times the pair's shared 1/<psil|psir> normalization factor, which
# cancels identically in the generalized Rayleigh quotient below, the
# only place this is used. inner() conjugates its first argument, so
# conj(wfl) is what turns it into the unconjugated wfl^T X wfr this
# convention calls for.
function nh_expect(wfl, X, wfr)
	return inner(conj(wfl)', X, wfr)
end


"""
    get_gs_generalized_nhdmrg(H, A, psi0; ...)

Non-Hermitian generalized eigenproblem: the generalized eigenvalue of
smallest real part, and the biorthogonal left/right eigenvector pair.
Returns `(lam, psil, psir)` with `<psil|psir> = 1` (ITensorNHDMRG's own
normalization, kept as-is).

Same outer loop as `get_gs_generalized` above with three changes, exactly
as in `pyitensor/nhdmrg.py`'s `nhdmrg_generalized()`: the inner solve is
ITensorNHDMRG's `nhdmrg` instead of `dmrg`, lambda is complex, and the
Rayleigh quotient is the *biorthogonal* one
`<psi_L|H|psi_R>/<psi_L|A|psi_R>` (evaluated through `nh_expect`, which
handles the left-vector convention -- see nhdmrg.jl's header). The
adjoint shift needs no separate bookkeeping here (unlike the pure-Python
port, which threads an explicit `HA` operator through): `nhdmrg` builds
the adjoint problem from its own `H` argument internally, and
`(H-lambda*A)^T = H^T - lambda*A` already because A is Hermitian *and*
(being a physical metric operator) symmetric, so the shift it applies to
its own transposed problem is the right one either way.

Each outer iteration warm-starts `nhdmrg` from the previous iteration's
*raw* (wfl,wfr) -- ITensorNHDMRG's own convention, which is what it
expects back -- and only converts to dmrgpy's convention for the
measurements and for the final return value.

`alg`/`biorthoalg` default to the same "onesided"+"fidelity" pair
`groundstate.py` uses for plain non-Hermitian ground states on this
backend.
"""
function get_gs_generalized_nhdmrg(H, A, psi0; nsweeps = 10, cutoff = 1e-8,
		maxm = 80, lam0 = NaN + 0.0im, alg = "onesided",
		biorthoalg = "fidelity", eigsolve_krylovdim = 30,
		eigsolve_maxite = 3)
	wfl = psi0
	wfr = psi0
	lam = ComplexF64(lam0)
	if isnan(real(lam))
		a0 = nh_expect(wfl, A, wfr)
		h0 = nh_expect(wfl, H, wfr)
		lam = abs(a0) > 1e-14 ? h0 / a0 : 0.0 + 0.0im
	end
	sweeps = make_sweeps(1, maxm, cutoff)
	for i = 1:nsweeps
		Heff = shifted_mpo(H, A, lam)
		# nhdmrg_solve, not a bare nhdmrg call: it carries nhdmrg.jl's
		# real-part-degeneracy tie-break, which matters just as much here
		# (Heff's spectrum is H's shifted by -lambda*A, so a conjugate
		# pair tying for the smallest real part stays a conjugate pair).
		# Its returned eigenvalue is unused -- lambda comes from the
		# Rayleigh quotient below, built from the *unrotated* H and A, so
		# the tie-break cannot leak into the result either way.
		_e, wfl, wfr = nhdmrg_solve(Heff, wfl, wfr, sweeps;
			alg = alg, biorthoalg = biorthoalg,
			eigsolve_krylovdim = eigsolve_krylovdim,
			eigsolve_maxite = eigsolve_maxite)
		a_psi = nh_expect(wfl, A, wfr)
		check_metric_norm(a_psi, "get_gs_generalized_nhdmrg")
		lam = nh_expect(wfl, H, wfr) / a_psi
	end
	psil, psir = nh_biorthogonal_pair(wfl, wfr)
	return lam, psil, psir
end
