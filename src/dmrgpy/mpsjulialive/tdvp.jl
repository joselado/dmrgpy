
# Real-time TDVP evolution, run natively in Julia (ITensorMPS.jl already
# exports tdvp() -- no separate ITensorTDVP.jl dependency needed). One
# call per full nt-step trajectory (not one Python<->Julia round trip per
# step), same design as kpm.jl's moment recursion. Mirrors
# pyitensor/chain.py's evolve_and_measure_tdvp()/quench_tdvp() (which in
# turn mirror mpscpp3/chain_session.h's C++ TDVP): one real-time step is
# exp(-i*dt*H), i.e. tdvp's evolution parameter is -im*dt.

function tdvp_step(H,psi,dt,cutoff,maxdim)
	# psi isn't always straight out of dmrg() or a previous tdvp_step()
	# call (e.g. evolution_ABA()/tdz.py's first step both start from an
	# operator applied via mpsjulialive/mpo.py's MPO.__mul__, which goes
	# through mpsalgebra.jl's applyoperator() -- confirmed directly to
	# leave a stale nonzero prime level on the result's Link indices).
	# tdvp()'s internal environment bookkeeping (ProjMPO) can't handle
	# that: it aborts deep inside KrylovKit's expintegrator ("... but the
	# indices are not permutations of each other", or a related
	# broadcast-shape failure in the same call chain). orthogonalize!()
	# alone does NOT clear the stale prime (confirmed directly), so the
	# Link indices need an explicit noprime() first; copy() so this never
	# mutates the caller's own MPS in place. Sanitizing on every step
	# rather than just the first is deliberate -- cheap next to the sweep
	# itself, and keeps every caller of tdvp_step correct by construction
	# instead of relying on each call site to remember to pre-clean.
	psi = noprime(copy(psi),"Link")
	orthogonalize!(psi,1)
	# No mindim override: an earlier version forced mindim=2 as a guessed
	# workaround for a DimensionMismatch that turned out to be the
	# stale-prime issue above, not a degenerate dim-1 bond -- confirmed
	# never actually necessary once that real cause was fixed. Neither
	# mpscpp3/chain_session.h's own tdvp_step (no MinDim set, so ITensor's
	# own default of 1) nor pyitensor/tdvp.py's svd calls (mindim=1
	# default) force a floor either; matching them keeps this backend
	# numerically comparable to both on small/near-product-state chains,
	# where forcing mindim=2 would otherwise pad in a spurious near-zero
	# singular value at a bond whose true Schmidt rank is 1.
	return tdvp(H,-im*dt,psi;time_step=-im*dt,cutoff=cutoff,maxdim=maxdim,
	            normalize=false,outputlevel=0)
end

function evolve_and_measure_tdvp(Hmpo,Aop,wf,nt,dt,cutoff,maxdim)
	# tdvp_step()'s internal normalize=false is deliberate and must stay
	# false at that shared level: advance_complex_time_step() (tdz.jl's
	# TDZ caller, see tdvp_step's own docstring) reuses tdvp_step for
	# complex-time evolution, where the norm is *meant* to decay -- that
	# decay is the whole damping mechanism TDZ's method relies on
	# (arXiv:2311.10909), not a numerical artifact to correct away.
	# Renormalizing there would silently break TDZ's physics. For plain
	# real-time evolution (this function), norm drift from per-step MPS
	# truncation *is* purely numerical, so renormalize explicitly here
	# instead -- same pattern quench_tdvp already uses below, and matches
	# mpscpp3/chain_session.h's own tdvp_step, which passes
	# "DoNormalize",true internally for exactly the real-time case.
	psi = wf*(1.0/sqrt(abs(inner(wf,wf))))
	correlator = ComplexF64[]
	for it=1:nt
		psi = tdvp_step(Hmpo,psi,dt,cutoff,maxdim)
		psi = psi*(1.0/sqrt(abs(inner(psi,psi))))
		push!(correlator,inner(psi,Aop,psi))
	end
	return correlator,psi
end

function apply_clean(A,psi0,maxdim,cutoff;alg="densitymatrix")
	# Deliberately NOT mpsalgebra.jl's applyoperator(), and deliberately
	# apply() rather than contract(): confirmed directly (linkinds(...)
	# checked explicitly) that contract(A,psi0) -- what applyoperator()
	# and its "iden_op fixes a bug" follow-up both build on -- leaves a
	# stale nonzero prime level on the result's Link indices, which
	# tdvp()'s internal environment bookkeeping (ProjMPO) can't handle:
	# it aborts deep inside KrylovKit's expintegrator with "... but the
	# indices are not permutations of each other" on the affected Link
	# index. apply() is ITensorMPS's own higher-level MPO-application
	# function (contract() is the lower-level primitive it and
	# applyoperator() are both built on) and does not have this problem.
	# orthogonalize!() gives the result a well-defined orthogonality
	# center, which tdvp() (like dmrg()) requires of its input MPS.
	psi1 = apply(A,psi0;maxdim=maxdim,cutoff=cutoff,alg=alg)
	orthogonalize!(psi1,1)
	return psi1
end

function quench_tdvp(Hshiftmpo,A1mpo,A2mpo,wf0,nt,dt,cutoff,maxdim)
	psi1 = apply_clean(A1mpo,wf0,maxdim,cutoff)
	psi2 = apply_clean(A2mpo,wf0,maxdim,cutoff)
	norm0 = sqrt(abs(inner(psi1,psi1)))
	correlator = ComplexF64[]
	for it=1:nt
		psi1 = tdvp_step(Hshiftmpo,psi1,dt,cutoff,maxdim)
		psi1 = psi1*(norm0/sqrt(abs(inner(psi1,psi1))))
		push!(correlator,inner(psi2,psi1))
	end
	return correlator,psi1
end
