# Site tensors of an MPS and of an MPO as dense arrays, for the
# tangent-space Krylov method (submode="TASK"), which runs on numpy arrays
# (pyitensor/task.py). Same layout as mpscpp3/bindings.cc's
# mps_site_arrays/mpo_site_arrays and pyitensor/chain.py's
# _mps_arrays_lpr/task.mpo_arrays: (left link, site, right link) for an
# MPS and (left link, primed site, site, right link) for an MPO, with
# length-1 axes at the open ends. As on v3 the tensors are exported as
# they are, with no change of gauge or norm. QN blocks are made dense.

function _dense_tensor(T)
	return hasqns(T) ? dense(T) : T
end

function _dense_index(i)
	return hasqns(i) ? removeqns(i) : i
end

function _site_array(T,inds_order)
	# inds_order: indices or `nothing` for a length-1 open end
	T = _dense_tensor(T)
	live = [_dense_index(i) for i in inds_order if i !== nothing]
	A = Array(T,live...)
	shape = [i === nothing ? 1 : dim(i) for i in inds_order]
	return ComplexF64.(reshape(A,shape...))
end

function mps_site_arrays(psi)
	n = length(psi)
	out = Array{ComplexF64,3}[]
	for i=1:n
		l = i > 1 ? linkind(psi,i-1) : nothing
		r = i < n ? linkind(psi,i) : nothing
		s = siteind(psi,i)
		push!(out,_site_array(psi[i],[l,s,r]))
	end
	return out
end

function mpo_site_arrays(H)
	n = length(H)
	out = Array{ComplexF64,4}[]
	for i=1:n
		l = i > 1 ? linkind(H,i-1) : nothing
		r = i < n ? linkind(H,i) : nothing
		s = firstind(H[i];tags="Site",plev=0)
		push!(out,_site_array(H[i],[l,prime(s),s,r]))
	end
	return out
end

function mps_mpo_mps(a,M,b)
	# <a|M|b>, a conjugated
	return inner(a',M,b)
end
