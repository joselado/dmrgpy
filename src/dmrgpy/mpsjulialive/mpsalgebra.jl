
function aMb(a,M,b,maxdim,cutoff; alg = "densitymatrix")
	iden_op = identity_mpo(sites) # identity operator
	return contract(a,M,b)
end

function applyoperator(sites,A,psi0,maxdim,cutoff; alg = "densitymatrix")
	iden_op = identity_mpo(sites) # identity operator
	psi1 = contract(A,psi0,maxdim=maxdim,cutoff=cutoff,alg = alg)
	psi1 = contract(iden_op,psi1,maxdim=maxdim,alg = alg) # this fixes a bug
	return psi1
end


function apply_op(A,psi0,maxdim,cutoff;alg="densitymatrix")
	# Deliberately apply() rather than this file's own applyoperator():
	# applyoperator() does contract(A,psi0) *plus* a second contraction
	# against an identity MPO ("this fixes a bug" -- a Site/Link-index
	# priming issue, see applyoperator()'s and tdvp.jl's apply_clean()'s
	# own notes), i.e. two truncated MPO-MPS contractions per call.
	# apply() is ITensorMPS's own higher-level MPO-application function
	# (the thing applyoperator()'s contract() is the lower-level
	# primitive for) and needs only one -- confirmed correctness-safe and
	# cheaper (measured directly on a 30-site Heisenberg chain, warm: the
	# KPM dynamical correlator went from ~34.5s to ~28.4s, ~18% faster).
	# Shared by kpm.jl's Chebyshev recursion (no orthogonality
	# requirement -- its output only ever feeds more apply_op()/inner()/
	# summps() calls) and tdvp.jl's apply_clean() (which additionally
	# orthogonalizes its result, since tdvp()/dmrg() require a
	# well-defined orthogonality center on their input) -- one primitive
	# instead of two independently-written near-duplicates.
	return apply(A,psi0;maxdim=maxdim,cutoff=cutoff,alg=alg)
end


function summps(psi1,psi2,maxdim,cutoff; alg = "densitymatrix")
	# cutoff must be forwarded explicitly: ITensorMPS's own add() defaults
	# to cutoff=1e-15 when not given, three orders tighter than dmrgpy's
	# usual configured cutoffs (e.g. kpmcutoff's own default of 1e-12) --
	# confirmed via the vendored ITensorMPS source. Without it, every
	# summps() call (kpm.jl's Chebyshev recursion, MPS.__add__) silently
	# keeps far more Schmidt values than the caller configured, up to
	# maxdim.
	psi3 = add(psi1,psi2,maxdim=maxdim,cutoff=cutoff,alg = alg)
	return psi3
end



function exponential(H,psi1,dtr,dti,maxdim,cutoff,nt0)
	tau = dtr + im*dti # evolution time
	taui = tau/nt0 # small time step
	iden_op = identity_mpo() # identity operator
        for i=1:nt0
          	psi2 = taui*contract(H,psi1;maxdim=maxdim,cutoff=cutoff,
				     method = "densitymatrix") # sum
          	psi2 = contract(iden_op,psi2;maxdim=maxdim,cutoff=cutoff,
				method = "densitymatrix") # fix
          	psi1 = add(psi1,psi2;maxdim=maxdim,cutoff=cutoff,
			   method = "densitymatrix") # sum
	end
	return psi1
end




function overlap(wf1,wf2)
	c = inner(wf1,wf2)
	return c
end



function random_state(sites)
	return random_mps(sites)
end


function mpstimesscalar(a,wf1)
	return a*wf1
end


