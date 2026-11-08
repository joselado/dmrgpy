# A^-1 b by BiCGSTAB over MPS, a line-by-line port of
# pyitensor/chain.py::_bicstab (and mpscpp3's bicstab). The inner products
# are the bilinear ones, inner(conj_mps(x),y) = sum_i x_i y_i, and the loop
# stops once sqrt(|Re r.r|) <= tol, with tol absolute: the caller
# (mpsalgebra.py::applyinverse_dmrg) passes delta*||b||. Every MPO
# application and every sum is truncated at maxdim and cutoff.
function mps_bicgstab(A,b,tol,max_it,maxdim,cutoff)
	x = copy(b)
	r_old = summps(b,(-1.0)*apply_op(A,x,maxdim,cutoff),maxdim,cutoff)
	r_ = copy(r_old)
	p = copy(r_old)
	k = 0
	while k < max_it
		Ap = apply_op(A,p,maxdim,cutoff)
		alpha = inner(conj_mps(r_old),r_)/inner(conj_mps(Ap),r_)
		s = summps(r_old,(-alpha)*Ap,maxdim,cutoff)
		As = apply_op(A,s,maxdim,cutoff)
		w = inner(conj_mps(As),s)/inner(conj_mps(As),As)
		x = summps(x,summps(alpha*p,w*s,maxdim,cutoff),maxdim,cutoff)
		r_new = summps(s,(-w)*As,maxdim,cutoff)
		res = sqrt(abs(real(inner(conj_mps(r_new),r_new))))
		if res <= tol
			break
		end
		beta = (alpha/w)*inner(conj_mps(r_new),r_)/inner(conj_mps(r_old),r_)
		p = summps(r_new,beta*summps(p,(-w)*Ap,maxdim,cutoff),maxdim,cutoff)
		r_old = r_new
		k += 1
	end
	return x
end
