
# Kernel Polynomial Method: Chebyshev-moment recursion run entirely in
# Julia (no per-step Python round trip), driven through mpsalgebra.jl's
# apply_op()/summps()/overlap(). Mirrors pyitensor/chain.py's own
# pure-Python _kpm_moments_full/_accelerated, and
# mpsjulialive/dynamics.py's earlier per-step Python orchestration of the
# same recursion (kept only as the Python-level fallback semantics --
# this file replaces the hot loop itself). apply_op() (rather than
# mpsalgebra.jl's applyoperator(), which does an extra, redundant MPO
# contraction per call -- see its own docstring) roughly halves the cost
# of this recursion's dominant work: measured directly on a 30-site
# Heisenberg chain, warm, the KPM dynamical correlator went from ~34.5s
# to ~28.4s (~18% faster), narrowing the gap to the compiled ITensor v3
# backend's ~23.4s from ~1.47x to ~1.21x slower.

function check_kpm_moment(lastmoment,bound)
	# Chebyshev moments of a correctly scaled Hamiltonian (spectrum
	# inside [-1,1]) satisfy |<vj|T_k|vi>| <= ||vi||*||vj|| = bound;
	# exponential growth beyond it means the scaled spectrum leaked
	# outside [-1,1] (band-edge estimate too tight for the chosen
	# kpm_scale) and every subsequent moment is garbage.
	if abs(lastmoment) > 1e3*(bound+1.0)
		error("KPM moments diverging: scaled spectrum outside [-1,1] "*
		      "(band-edge estimate too tight; increasing kpm_scale "*
		      "widens the safety margin)")
	end
end


function kpm_moments_full(jlmpo,vi,vj,n,kpmmaxm,kpmcutoff)
	v = 1.0*vi
	am = 1.0*vi
	a = apply_op(jlmpo,v,kpmmaxm,kpmcutoff)
	# legitimate moments <vj|T_k|vi> are bounded by ||vi||*||vj|| (NOT by
	# the zeroth moment <vj|vi>, which can be ~0 for a near-orthogonal
	# cross-correlator pair)
	bound = sqrt(abs(real(inner(vi,vi))*real(inner(vj,vj))))
	out = ComplexF64[]
	push!(out,inner(vj,v))
	push!(out,inner(vj,a))
	for k=1:n
		ap = apply_op(jlmpo,a,kpmmaxm,kpmcutoff)
		ap = summps(2.0*ap,(-1.0)*am,kpmmaxm,kpmcutoff)
		push!(out,inner(vj,ap))
		check_kpm_moment(out[end],bound)
		am = 1.0*a
		a = 1.0*ap
	end
	return out
end


function kpm_moments_accelerated(jlmpo,vi,n,kpmmaxm,kpmcutoff)
	a = apply_op(jlmpo,vi,kpmmaxm,kpmcutoff)
	am = 1.0*vi
	mu0 = inner(vi,vi)
	mu1 = inner(vi,a)
	bound = abs(mu0) # here vi==vj, so the moment bound ||vi||*||vj|| is just mu0
	out = ComplexF64[mu0,mu1]
	for k=1:div(n,2)
		ap = apply_op(jlmpo,a,kpmmaxm,kpmcutoff)
		ap = summps(2.0*ap,(-1.0)*am,kpmmaxm,kpmcutoff)
		bk = 2.0*inner(a,a)-mu0
		bk1 = 2.0*inner(a,ap)-mu1
		push!(out,bk)
		push!(out,bk1)
		check_kpm_moment(out[end],bound)
		am = 1.0*a
		a = 1.0*ap
	end
	return out
end


function same_mps(vi,vj,maxm,cutoff)
	d = summps(1.0*vi,(-1.0)*vj,maxm,cutoff)
	dd = sqrt(abs(real(inner(d,d))))
	return dd<1e-10
end
