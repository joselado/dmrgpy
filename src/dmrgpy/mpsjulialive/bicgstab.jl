# A^-1 b by BiCGSTAB over MPS, a line-by-line port of
# pyitensor/chain.py::_bicstab (and mpscpp3's bicstab). The inner products
# are the bilinear ones, inner(conj_mps(x),y) = sum_i x_i y_i, and the loop
# stops once sqrt(|Re r.r|) <= tol, with tol absolute: the caller
# (mpsalgebra.py::applyinverse_dmrg) passes delta*||b||. Every MPO
# application and every sum is truncated at maxdim and cutoff.
#
# The sums are a direct sum compressed by SVD (bicg_sum), not summps,
# whose density-matrix algorithm squares the weights and so resolves a
# sum only to about 1e-8 of its terms. The residual updates here are
# exactly such cancellations, and with summps the recursion broke down:
# CVM_explicit of (Sx0,Sy1) on the 6-site Heisenberg chain, exactly 0,
# came back as 4.8e9 at omega=1.8 at the default cutoff 1e-12.
function bicg_sum(a,b,maxdim,cutoff)
	return truncate(add(a,b;alg="directsum");maxdim=maxdim,cutoff=cutoff)
end

function mps_bicgstab(A,b,tol,max_it,maxdim,cutoff)
	x = copy(b)
	r_old = bicg_sum(b,(-1.0)*apply_op(A,x,maxdim,cutoff),maxdim,cutoff)
	r_ = copy(r_old)
	p = copy(r_old)
	k = 0
	while k < max_it
		Ap = apply_op(A,p,maxdim,cutoff)
		alpha = inner(conj_mps(r_old),r_)/inner(conj_mps(Ap),r_)
		s = bicg_sum(r_old,(-alpha)*Ap,maxdim,cutoff)
		As = apply_op(A,s,maxdim,cutoff)
		w = inner(conj_mps(As),s)/inner(conj_mps(As),As)
		x = bicg_sum(x,bicg_sum(alpha*p,w*s,maxdim,cutoff),maxdim,cutoff)
		r_new = bicg_sum(s,(-w)*As,maxdim,cutoff)
		res = sqrt(abs(real(inner(conj_mps(r_new),r_new))))
		if res <= tol
			break
		end
		beta = (alpha/w)*inner(conj_mps(r_new),r_)/inner(conj_mps(r_old),r_)
		p = bicg_sum(r_new,beta*bicg_sum(p,(-w)*Ap,maxdim,cutoff),maxdim,cutoff)
		r_old = r_new
		k += 1
	end
	return x
end
