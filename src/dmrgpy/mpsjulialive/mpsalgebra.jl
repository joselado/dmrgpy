
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


function summps(psi1,psi2,maxdim; alg = "densitymatrix")
	psi3 = add(psi1,psi2,maxdim=maxdim,alg = alg)
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


