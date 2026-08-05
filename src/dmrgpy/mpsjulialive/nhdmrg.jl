
# Non-Hermitian DMRG on the live Julia backend, via the real
# ITensorNHDMRG.jl package (declared in ../juliapkg.json) rather than one
# of this codebase's own ports of it. get_gs.jl's get_gs_nhdmrg already
# calls it for plain non-Hermitian ground-state energies and keeps only
# the right eigenvector; this file exposes the full biorthogonal pair, in
# the convention the rest of dmrgpy uses, and handles the two places
# where "the rest of dmrgpy" and ITensorNHDMRG genuinely differ.
#
# (1) LEFT-VECTOR CONVENTION. ITensorNHDMRG's "adjoint" problem is the
# TRANSPOSE one, not the conjugate-transpose one: its ProjNHMPO
# constructor builds the second projector from
# `dag(swapprime(conj(H), 0 => 1))`, and since dag() is itself a
# conjugation here (these sites carry no QNs, see sites.jl's plain
# Index(...) construction, so dag() has no arrows to flip), the two
# conjugations cancel and what it actually sweeps against is
# swapprime(H, 0 => 1) == H^T. So the returned `wfl` satisfies
#
#     H^T |wfl> = lambda |wfl>,   with inner(wfl,wfr) == 1
#
# whereas every other dmrgpy backend (nhdmrg.py's documented convention,
# and what its residual certificate checks) defines the left eigenvector
# by
#
#     H^dagger |psil> = conj(lambda) |psil>,  with <psil|psir> == 1.
#
# The two differ by exactly a complex conjugation: conj() of the first
# equation is H^dagger conj(wfl) = conj(lambda) conj(wfl). Hence
# psil = conj(wfl), renormalized because inner() conjugates its first
# argument and so the ITensorNHDMRG normalization does not survive the
# conjugation.
#
# This is NOT visible on a complex-*symmetric* Hamiltonian (H^T == H),
# where the transpose and adjoint left vectors coincide up to that same
# conjugation -- which covers both of tests/test_nh_dmrg.py's original
# models (the staggered-imaginary-potential fermionic chain and the
# PT-symmetric XX chain are both complex symmetric), so a genuinely
# non-symmetric model is what pins this down (see tests/test_nh_dmrg.py's
# nh_asymmetric_hopping_chain, added for exactly this reason). Confirmed
# directly: without the conjugation the left residual
# ||H^dagger|psil> - conj(lambda)|psil>|| sits at ~2 while the right
# residual is ~1e-15.
#
# (2) UNANCHORED LEFT SOLVE ON A REAL-PART-DEGENERATE SPECTRUM. dmrgpy
# targets the eigenvalue of smallest *real* part. When two eigenvalues
# tie for it -- a complex-conjugate pair a±ib, which is the generic
# situation in PT-symmetric and Hatano-Nelson-like models -- nothing in
# ITensorNHDMRG ties its left solve to whichever member its right solve
# picked, so the two can converge to *different* eigenvalues. dmrgpy's
# own ports anchor against exactly this (mpscpp3/chain_session.h's
# arnoldi_smallest_real and its Sel comment). Observed directly, and
# deterministically (not a bad random draw -- every attempt from every
# start did it) on an asymmetric-hopping chain: conj(wfl) came back
# satisfying H^dagger conj(wfl) = lambda conj(wfl) instead of
# conj(lambda) conj(wfl), i.e. it was the left vector of the conjugate
# partner, and its overlap with wfr was exactly 0 (the two belong to
# different eigenvalues).
#
# nhdmrg_solve() below breaks the tie instead of anchoring, which needs
# no change inside ITensorNHDMRG: solving for exp(i*theta)*H leaves every
# eigenVECTOR untouched and maps every eigenVALUE to
# exp(i*theta)*lambda, so Re(exp(i*theta)*(a±ib)) = a*cos(theta) -+
# b*sin(theta) is no longer degenerate, and the eigenvalue is mapped back
# exactly by dividing by the same phase.
#
# The rotation is NOT harmless on its own, and an earlier version of this
# file wrongly claimed it was ("a theta this small cannot change which
# eigenvalue is selected"). It shifts every eigenvalue's real part by
# -Im(lambda)*sin(theta), which grows with |Im lambda| and therefore with
# system size and non-Hermitian coupling strength: at n=4 with gamma=0.7
# the imaginary parts are ~0.5 and the shift is ~5e-3, but at n=20 they
# reach ~7 and the shift is ~0.07 -- comparable to a real low-lying
# real-part gap. A large enough theta can therefore make ITensorNHDMRG
# converge on a genuinely DIFFERENT eigenvalue, which then maps back to a
# perfectly valid eigenpair of H that is simply not the smallest-real-part
# one. Neither residual in nhdmrg.py's certificate can catch that (both
# sit at ~1e-15 for any true eigenpair), so it would be a silently wrong
# ground-state energy.
#
# Two things make the tie-break safe regardless of theta:
#
#   - the untied solve's own eigenvalue e0 is kept as the anchor. Only its
#     LEFT vector is suspect -- its right vector and eigenvalue are a
#     converged eigenpair (verified downstream at ~1e-15) -- so a
#     tie-break result is accepted only if it reproduces e0. Anything else
#     means the rotation moved the target, and is rejected rather than
#     returned.
#   - theta is tried LARGEST-first, which is not the obvious order and was
#     established by measurement. Breaking the tie needs only theta>0 in
#     exact arithmetic, but numerically the induced real-part split has to
#     be large enough for the solver to actually resolve, and below that
#     the run comes back nearly-degenerate and badly converged rather than
#     wrong: starting the ladder at 1e-4 reproducibly returned a pair with
#     the right eigenvalue but a ~3e-2 residual on a chain where 1e-2 gives
#     ~1e-15. So 1e-2 leads, and the smaller angles are fallbacks for the
#     case the anchor check rejects it (large |Im lambda|, where a big
#     rotation is what risks re-targeting in the first place).
#
# The two conditions fail in opposite, and safe, directions: too small a
# theta gives a poorly-converged pair, which nhdmrg.py's residual
# certificate catches and retries; too large a theta gives a
# well-converged pair for the wrong eigenvalue, which only the anchor
# check can catch. That asymmetry is why the anchor check is mandatory and
# the angle choice is merely tuning.
#
# If no angle satisfies both conditions, this gives up and lets
# nh_biorthogonal_pair raise, which nhdmrg.py treats as a failed attempt.
# Refusing to answer is the right outcome there: the alternative is
# returning the wrong eigenvalue with a clean bill of health.

using ITensors
using ITensorMPS
using ITensorNHDMRG


# The bilinear (unconjugated) overlap wfl^T*wfr, which is what
# <psil|psir> becomes under psil = conj(wfl) since inner() conjugates its
# first argument. Zero exactly when the two vectors belong to different
# eigenvalues -- the signature of the unanchored-left-solve problem
# above, and the quantity the biorthogonal normalization divides by.
function nh_pair_overlap(wfl, wfr)
	return inner(conj(wfl), wfr)
end


# The single point of contact with ITensorNHDMRG's own solver, shared by
# get_gs.jl's get_gs_nhdmrg (plain non-Hermitian ground-state energies)
# and by nhdmrg_solve below -- so its call signature and defaults exist in
# exactly one place. Two independent copies used to drift apart silently
# whenever ITensorNHDMRG changed a keyword or a default, leaving
# chain.gs_energy() and chain.nhdmrg() running differently-configured
# solvers on the same chain.
function nhdmrg_raw(H, psil0, psir0, sweeps; alg = "onesided",
		biorthoalg = "biorthoblock", eigsolve_krylovdim = 30,
		eigsolve_maxite = 3)
	return run_quiet() do
		nhdmrg(
			H,
			psil0,
			psir0,
			sweeps;
			alg = alg,
			biorthoalg = biorthoalg,
			outputlevel = 1,
			eigsolve_krylovdim = eigsolve_krylovdim,
			eigsolve_maxiter = eigsolve_maxite,
		)
	end
end


# Rotation angles tried by the tie-break, largest first (see the header --
# the order is measured, not obvious).
const TIEBREAK_ANGLES = [1.0e-2, 1.0e-3, 1.0e-4]

# 1e-8, not ~1e-14: the collapse is total when it happens (measured at
# ~1e-15 against a healthy value of ~0.5), so anything in between is
# comfortably a threshold rather than a tuned constant.
const PAIR_OVERLAP_TOL = 1.0e-8

# How far a tie-break run's eigenvalue may sit from the untied run's and
# still count as "the same eigenvalue". Has to straddle a wide gap: two
# runs that genuinely found the same eigenvalue have differed by up to
# ~1e-5 here (ordinary DMRG re-convergence noise, measured), while two
# *different* eigenvalues of these models differ by O(0.1-1). 1e-3
# relative sits comfortably between the two, so it is a separator rather
# than a tuned constant.
same_eigenvalue(e, e0) = abs(e - e0) <= max(1.0e-6, 1.0e-3 * (1 + abs(e0)))


"""
    nhdmrg_solve(H, psil0, psir0, sweeps; tiebreak, ...)

One ITensorNHDMRG sweep schedule, returning its raw `(lambda,wfl,wfr)` in
ITensorNHDMRG's own convention, but with the real-part-degeneracy
tie-break of this file's header applied when needed: if the left and
right vectors come back belonging to different eigenvalues (bilinear
overlap ~0), the solve is retried against `exp(i*theta)*H` over
`TIEBREAK_ANGLES`, and the first result that both fixes the overlap AND
reproduces the untied run's eigenvalue is returned.

If no angle achieves both, the untied result is returned unchanged, whose
collapsed overlap makes `nh_biorthogonal_pair` raise -- deliberately, see
the header: returning a tie-break result that landed on a different
eigenvalue would be a silently wrong answer, and no downstream check
could catch it. `tiebreak = false` disables the whole mechanism.
"""
function nhdmrg_solve(H, psil0, psir0, sweeps; tiebreak = true, kwargs...)
	e0, wfl, wfr = nhdmrg_raw(H, psil0, psir0, sweeps; kwargs...)
	if !tiebreak || abs(nh_pair_overlap(wfl, wfr)) > PAIR_OVERLAP_TOL
		return e0, wfl, wfr
	end
	for theta in TIEBREAK_ANGLES
		phase = exp(im * theta)
		e, l, r = nhdmrg_raw(phase * H, psil0, psir0, sweeps; kwargs...)
		e = e / phase # exact: the rotation scales every eigenvalue by phase
		if abs(nh_pair_overlap(l, r)) > PAIR_OVERLAP_TOL && same_eigenvalue(e, e0)
			return e, l, r
		end
	end
	return e0, wfl, wfr
end


"""
    nh_biorthogonal_pair(wfl, wfr)

Convert ITensorNHDMRG's (wfl,wfr) output into dmrgpy's left/right
convention: returns `(psil, psir)` with `H^dagger|psil> =
conj(lambda)|psil>` and `<psil|psir> = 1`. See this file's header for why
the conjugation is needed.
"""
function nh_biorthogonal_pair(wfl, wfr)
	psil = conj(wfl)
	s = nh_pair_overlap(wfl, wfr) # == <psil|psir>, before renormalization
	if !(abs(s) > 1e-14)
		# Reachable only if nhdmrg_solve's tie-break didn't help either
		# (or was disabled); dividing through would hand back an Inf/NaN
		# MPS that looks like an ordinary result. nhdmrg.py's driver
		# treats this as a failed attempt and redraws.
		error("nh_biorthogonal_pair: <psil|psir> collapsed to ~0 -- the " *
		      "left/right pair is not biorthogonalizable")
	end
	# inner(c*psil, psir) == conj(c)*inner(psil,psir), so the scale factor
	# that lands <psil|psir> exactly on 1 is 1/conj(s), not 1/s.
	return psil / conj(s), wfr
end


"""
    get_nhdmrg_pair(H, psi0; ...)

One NH-DMRG run from `psi0`, returning `(lambda, psil, psir)` in dmrgpy's
convention (see `nh_biorthogonal_pair`). Same solver and same defaults as
get_gs.jl's `get_gs_nhdmrg`, which stays as it is -- it feeds
groundstate.py's plain non-Hermitian `gs_energy()` path, which only ever
wanted the right eigenvector.
"""
# Deliberately NOT given self.noise, unlike the Hermitian
# get_gs_generalized: DMRG's density-matrix noise term is defined against
# a Hermitian density matrix, and ITensorNHDMRG's biorthogonal truncation
# (the "fidelity" algorithm's rho=(rho_l+rho_r)/2 isometry) is not that.
# Feeding a noisy Sweeps through it was tried and measurably broke
# previously-converged runs -- the tie-break's own anchor check started
# rejecting every angle on a chain that had converged to ~1e-15 without
# it. The session backends do pass self.noise down their NH-DMRG path,
# so this is a real backend difference; it is documented rather than
# forced, because forcing it makes this backend worse, not more faithful.
function get_nhdmrg_pair(H, psi0; nsweeps = 10, cutoff = 1e-8, maxm = 80,
		alg = "onesided", biorthoalg = "biorthoblock",
		eigsolve_krylovdim = 30, eigsolve_maxite = 3, tiebreak = true)
	sweeps = make_sweeps(nsweeps, maxm, cutoff)
	e, wfl, wfr = nhdmrg_solve(H, psi0, psi0, sweeps; tiebreak = tiebreak,
		alg = alg, biorthoalg = biorthoalg,
		eigsolve_krylovdim = eigsolve_krylovdim,
		eigsolve_maxite = eigsolve_maxite)
	psil, psir = nh_biorthogonal_pair(wfl, wfr)
	return e, psil, psir
end
