#include("common.jl")
#using ITensors
#include("read_operator.jl")
#include("read_wf.jl")
#include("get_input.jl")
#include("get_sites.jl")
#include("get_gs.jl")
#include("write_in_file.jl")

function dynamical_correlator_kpm()
    sites = get_sites()
    n,wf,H = get_num_polynomials_and_hamiltonian(sites) # get several objects
    A = MPO(read_operator("kpm_multioperator_i.in"),sites) # operator A
    B = MPO(read_operator("kpm_multioperator_j.in"),sites) # operator B
    get_moments(A,B,H,wf,n) # compute the moments
end

function get_num_polynomials_and_hamiltonian(sites)
	"""Compute bandwidth of the Hamiltonian"""
	ampo = read_operator("hamiltonian.in")
	H = MPO(ampo,sites) # Hamiltonian
	Hm = -1*H # minus Hamiltonian
	psi0 = randomMPS(sites)
	psi1 = randomMPS(sites)
        sweeps = Sweeps(get_int("nsweeps"))
        maxdim!(sweeps, get_int("maxm"))
        cutoff!(sweeps, get_float("cutoff"))
	e0, wf = dmrg(H,psi0, sweeps) # lower energy
        e1, psi = dmrg(Hm,psi1, sweeps) # upper energy
	e1 = -e1 # upper energy
	n = (e1-e0)/get_float("kpm_delta")
	n = n*get_int("kpm_n_scale") # scaling of the KPM method
	n = trunc(Int, n) # return the number of polynomials
	write_in_file("KPM_NUM_POLYNOMIALS.OUT",n) # number of polynomials
	ampo += (-(e1+e0)/2.0,"Id",1) # shift the Hamiltonian
	Ho = MPO(ampo,sites) # create shifted Hamiltonian
	scale = (e1-e0)*get_float("kpm_scale") # scale for the hamiltonian
	scale = 1.0/scale # inverse scale
	write_in_file("KPM_SCALE.OUT",[e0,e1,scale])
	write_in_file("GS_ENERGY.OUT",e0)
	return n,wf,scale*Ho # return scaled Hamiltonian
end

function get_scaled_hamiltonian(sites)
	"""Get the reescaled Hamiltonian"""
	ampo = read_operator("hamiltonian.in")
	H = MPO(ampo,sites) # Hamiltonian
end

function get_moments(A,B,H,wf,n)
	"""Compute the moments using the KPM recursion"""
	maxdim = get_int("maxm")
	vi = applyMPO(A,wf) # first vector
	vj = applyMPO(B,wf) # second vector
	a = applyMPO(H,vi) # first vector
	am = 1.0*vi # define
	bk = inner(vj,vi) # scalar product
	bk1 = inner(vj,a) # scalar product
	write_in_file("KPM_MOMENTS.OUT",bk,"w")
	write_in_file("KPM_MOMENTS.OUT",bk1,"a")
	for i=1:n # loop over polynomials
	  @time bk,am,a = kpm_iterate(am,a,H,vj) # perform one iteration
#	  write_in_file("KPM_MOMENTS.OUT",bk,"a")
#	  truncate!(a,maxdim=maxdim) 
#	  truncate!(am,maxdim=maxdim) 
#	  write_in_file("KPM_MOMENTS.OUT",bk,"a")
        end
end

function kpm_iterate(am,a,H,vj)
	  """Perform a single iteration"""
	  ap = applyMPO(H,a) # first vector
	  ap = sum(2*ap,-1*am) # redefine the vector
	  bk = inner(vj,ap) # scalar product
#	  truncate!(a,maxdim=maxdim) 
#	  truncate!(am,maxdim=maxdim) 
	  return bk,a,ap # redefine
end

#dynamical_correlator_kpm()
