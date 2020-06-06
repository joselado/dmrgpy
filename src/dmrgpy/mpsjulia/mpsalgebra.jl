function applyoperator()
	sites = get_sites()
	nameop = get_string("applyoperator_multioperator")
	A = read_mpo(nameop) # read the operator
	psi0 = load_mps(get_string("applyoperator_wf0"))
	psi1 = contract(A,psi0,maxdim=get_int("maxm"))
	save_mps(get_string("applyoperator_wf1"),psi1)
end


function exponential()
	dtr = get_float("tevol_dt_real")
	dti = get_float("tevol_dt_imag")
	maxdim = get_int("maxm")
	H = read_mpo("hamiltonian.in")  # input wavefunction
	psi1 = load_mps("input_wavefunction.mps") # Operator to exponentiate
	nt0 = get_int("tevol_n") # number of steps
	tau = dtr + im*dti # evolution time
	taui = tau/nt0 # small time step
	iden_op = identity_mpo() # identity operator
        for i=1:nt0
          	psi2 = taui*contract(H,1*psi1,maxdim=maxdim) # sum
          	psi2 = contract(iden_op,psi2,maxdim=maxdim) # this fixes a bug
          	psi1 = add(1*psi1,1*psi2,maxdim=maxdim) # sum
	end
	save_mps("output_wavefunction.mps",psi1) # write result
end




function overlap()
	wf1 = load_mps("overlap_wf1.mps")
	wf2 = load_mps("overlap_wf2.mps")
	c = inner(wf1,wf2)
	write_in_file("OVERLAP.OUT",convert(Complex,c),"w")
end

