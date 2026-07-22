
# Real-time TDVP evolution, run natively in Julia (ITensorMPS.jl already
# exports tdvp() -- no separate ITensorTDVP.jl dependency needed). One
# call per full nt-step trajectory (not one Python<->Julia round trip per
# step), same design as kpm.jl's moment recursion. Mirrors
# pyitensor/chain.py's evolve_and_measure_tdvp()/quench_tdvp() (which in
# turn mirror mpscpp3/chain_session.h's C++ TDVP): one real-time step is
# exp(-i*dt*H), i.e. tdvp's evolution parameter is -im*dt.

function tdvp_step(H,psi,dt,cutoff,maxdim)
	# mindim=2 is a cheap defensive floor against degenerate dim-1 bonds;
	# kept even though the actual bug hit here (see apply_clean() below)
	# turned out to be elsewhere.
	return tdvp(H,-im*dt,psi;time_step=-im*dt,cutoff=cutoff,maxdim=maxdim,
	            mindim=2,normalize=false,outputlevel=0)
end

function evolve_and_measure_tdvp(Hmpo,Aop,wf,nt,dt,cutoff,maxdim)
	# wf isn't always straight out of dmrg() (e.g. evolution_ABA() first
	# applies an operator via mpsjulialive/mpo.py's MPO.__mul__, which
	# goes through mpsalgebra.jl's applyoperator() -- confirmed directly
	# to leave a stale nonzero prime level on the result's Link indices,
	# same as apply_clean()'s docstring below explains for quench_tdvp's
	# inputs. orthogonalize!() alone does NOT clear this -- confirmed
	# directly, the stale prime survives it -- so the Link indices need
	# an explicit noprime() first. copy() first so this doesn't mutate
	# the caller's own wf in place.
	psi = noprime(copy(wf),"Link")
	orthogonalize!(psi,1)
	correlator = ComplexF64[]
	for it=1:nt
		psi = tdvp_step(Hmpo,psi,dt,cutoff,maxdim)
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
