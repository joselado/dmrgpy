
# Excited states via ITensorMPS.jl's own orthogonality-penalty dmrg()
# method (dmrg(H::MPO, Ms::Vector{MPS}, psi0::MPS, sweeps; weight)) --
# mirrors pyitensor/chain.py's excited_states()/dmrg_excited(), which in
# turn mirrors mpscpp3/chain_session.h's Chain::excited_states().

function excited_state_dmrg(H,wfs,psi0,weight,nsweeps,cutoff,maxm)
	# make_sweeps/run_quiet are defined in get_gs.jl, loaded earlier by
	# juliasession.py's initialize() -- shared rather than a third
	# copy-pasted Sweeps-construction+stdout-silencing block.
	sweeps = make_sweeps(nsweeps,maxm,cutoff)
	run_quiet() do
		dmrg(H,wfs,psi0,sweeps;weight=weight)
	end
end


function excited_states_dmrg(H,psi0_gs,n,weight,sites,nsweeps,cutoff,maxm)
	# whole n-state loop runs in one Julia call, same design as
	# kpm.jl/tdvp.jl -- one new random-start warm state per additional
	# excited state, deflated against every wavefunction found so far
	# via dmrg()'s own orthogonality-penalty weight
	wfs = MPS[psi0_gs]
	for i=2:n
		psi0 = random_state(sites)
		e1,psi1 = excited_state_dmrg(H,wfs,psi0,weight,nsweeps,cutoff,maxm)
		push!(wfs,psi1)
	end
	# re-evaluate exact <wf|H|wf> for every state (including the ground
	# state) rather than reusing dmrg()'s own penalized eigenvalue
	# estimate -- mirrors mpscpp3/chain_session.h's Chain::excited_states
	energies = ComplexF64[inner(w,H,w) for w in wfs]
	return energies,wfs
end
