

function applyoperator(sites,A,psi0,maxdim,cutoff)
	iden_op = identity_mpo(sites) # identity operator
	psi1 = contract(A,psi0,maxdim=maxdim,cutoff=cutoff,alg = "densitymatrix")
	psi1 = contract(iden_op,psi1,maxdim=maxdim,alg = "densitymatrix") # this fixes a bug
	return psi1
end


function summps(psi1,psi2,maxdim)
	psi3 = add(psi1,psi2,maxdim=maxdim,alg = "densitymatrix")
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


