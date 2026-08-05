
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
# exactly by multiplying by exp(-i*theta). For a spectrum whose smallest
# real part is already unique and separated by a gap, a theta this small
# cannot change which eigenvalue is selected, so this is a pure tie-break
# rather than a change of target -- confirmed on the model above, where
# the returned eigenvalue is bit-for-bit what the untied solve returned
# and only the left vector changes (to the correct partner, left residual
# ~2 -> ~1e-15).

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


"""
    nhdmrg_solve(H, psil0, psir0, sweeps; tiebreak, ...)

One ITensorNHDMRG sweep schedule, returning its raw `(lambda,wfl,wfr)` in
ITensorNHDMRG's own convention, but with the real-part-degeneracy
tie-break of this file's header applied when needed: if the left and
right vectors come back belonging to different eigenvalues (bilinear
overlap ~0), the solve is repeated once against `exp(i*tiebreak)*H` and
the eigenvalue mapped back. `tiebreak = 0.0` disables it.
"""
function nhdmrg_solve(H, psil0, psir0, sweeps; tiebreak = 1e-2, kwargs...)
	e, wfl, wfr = nhdmrg_raw(H, psil0, psir0, sweeps; kwargs...)
	# 1e-8, not ~1e-14: the collapse is total when it happens (measured
	# at ~1e-15 against a healthy value of ~0.5), so anything in between
	# is comfortably a threshold rather than a tuned constant.
	if tiebreak == 0.0 || abs(nh_pair_overlap(wfl, wfr)) > 1e-8
		return e, wfl, wfr
	end
	phase = exp(im * tiebreak)
	e, wfl, wfr = nhdmrg_raw(phase * H, psil0, psir0, sweeps; kwargs...)
	return e / phase, wfl, wfr
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
function get_nhdmrg_pair(H, psi0; nsweeps = 10, cutoff = 1e-8, maxm = 80,
		alg = "onesided", biorthoalg = "biorthoblock",
		eigsolve_krylovdim = 30, eigsolve_maxite = 3, tiebreak = 1e-2)
	sweeps = make_sweeps(nsweeps, maxm, cutoff)
	e, wfl, wfr = nhdmrg_solve(H, psi0, psi0, sweeps; tiebreak = tiebreak,
		alg = alg, biorthoalg = biorthoalg,
		eigsolve_krylovdim = eigsolve_krylovdim,
		eigsolve_maxite = eigsolve_maxite)
	psil, psir = nh_biorthogonal_pair(wfl, wfr)
	return e, psil, psir
end
