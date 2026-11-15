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


# The inverse of mps_site_arrays: an MPS on the site indices `sites` from
# arrays of shape (left link, site, right link), length-1 axes at the open
# ends, as pyitensor/rsi.py returns them (mpsalgebratk/hadamard.py). New
# link indices are made for every bond. Dense sites only, as the export.
function mps_from_site_arrays(arrays,sites)
	n = length(sites)
	length(arrays) == n || error("mps_from_site_arrays: $(length(arrays)) "*
		"arrays for $n sites")
	links = [Index(size(arrays[i],3),"Link,l=$i") for i=1:n-1]
	psi = MPS(n)
	for i=1:n
		A = ComplexF64.(Array(arrays[i]))
		size(A,2) == dim(sites[i]) || error("mps_from_site_arrays: site "*
			"$i has dimension $(dim(sites[i])), the array $(size(A,2))")
		inds_i = Index[]
		shape = Int[]
		if i > 1
			size(A,1) == dim(links[i-1]) || error("mps_from_site_arrays: "*
				"bond $(i-1) does not match")
			push!(inds_i,links[i-1]); push!(shape,size(A,1))
		end
		push!(inds_i,sites[i]); push!(shape,size(A,2))
		if i < n
			push!(inds_i,links[i]); push!(shape,size(A,3))
		end
		psi[i] = ITensor(reshape(A,shape...),inds_i...)
	end
	return psi
end


# psi compressed by SVD to at most maxdim states per bond and the
# truncation cutoff, for the exact Kronecker product of
# mpsalgebratk/hadamard.py's method="direct"
function compress_mps(psi,maxdim,cutoff)
	return truncate(psi;maxdim=maxdim,cutoff=cutoff)
end
