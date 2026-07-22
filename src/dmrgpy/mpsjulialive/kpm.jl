
# Kernel Polynomial Method: Chebyshev-moment recursion run entirely in
# Julia (no per-step Python round trip), driven through apply_op() below
# plus summps/overlap from mpsalgebra.jl. Mirrors pyitensor/chain.py's
# own pure-Python _kpm_moments_full/_accelerated, and
# mpsjulialive/dynamics.py's earlier per-step Python orchestration of the
# same recursion (kept only as the Python-level fallback semantics --
# this file replaces the hot loop itself).

function apply_op(A,psi0,maxdim,cutoff;alg="densitymatrix")
	# Deliberately apply() rather than mpsalgebra.jl's applyoperator():
	# applyoperator() does contract(A,psi0) *plus* a second contraction
	# against an identity MPO ("this fixes a bug" -- a Site-index priming
	# issue), i.e. two truncated MPO-MPS contractions per call. apply()
	# is ITensorMPS's own higher-level MPO-application function (the
	# thing applyoperator()'s contract() is the lower-level primitive
	# for) and needs only one -- confirmed correctness-safe and cheaper
	# during the TDVP work (see tdvp.jl's apply_clean(), which switched
	# for the same reason to fix a different bug). The Chebyshev
	# recursion below calls this once per moment, so this cuts one of
	# the two contractions applyoperator() did per step -- though not a
	# clean 2x, since applyoperator()'s second contraction is against a
	# trivial bond-dimension-1 identity MPO, much cheaper than the first
	# contraction against the real (bond-growing) scaled Hamiltonian.
	# Measured directly on a 30-site Heisenberg chain, warm: the KPM
	# dynamical correlator went from ~34.5s to ~28.4s (~18% faster),
	# narrowing the gap to the compiled ITensor v3 backend's ~23.4s from
	# ~1.47x to ~1.21x slower.
	return apply(A,psi0;maxdim=maxdim,cutoff=cutoff,alg=alg)
end


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
		ap = summps(2.0*ap,(-1.0)*am,kpmmaxm)
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
		ap = summps(2.0*ap,(-1.0)*am,kpmmaxm)
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


function same_mps(vi,vj,maxm)
	d = summps(1.0*vi,(-1.0)*vj,maxm)
	dd = sqrt(abs(real(inner(d,d))))
	return dd<1e-10
end
