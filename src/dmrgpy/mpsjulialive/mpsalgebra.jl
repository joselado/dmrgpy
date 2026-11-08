
# <a|M|b>, contracted exactly: no MPO-MPS product is formed and nothing
# is truncated. inner() conjugates its first argument. The primed bra is
# the form ITensorMPS asks for, with M carrying primes 0 (in) and 1 (out).
function aMb(a,M,b)
	return inner(a',M,b)
end


# The complex conjugate of every site tensor. With dense indices dag is
# conj; with QN indices it also flips the arrows, which is what a bra
# needs. inner(conj_mps(x),y) is then the bilinear sum_i x_i y_i, the
# product pyitensor's and mpscpp3's BiCGSTAB use.
function conj_mps(psi)
	return dag(psi)
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



# exp(tau*H)|psi>, as nt0 steps of the second-order Taylor operator
# E = 1 + z H + z^2/2 H^2 with z = tau/nt0. The same construction as
# pyitensor/chain.py::_custom_exp and mpscpp3's custom_exp, so that the
# three backends agree at equal nt0. E is built exactly (direct sums) and
# then compressed at a relative cutoff of 1e-28, which drops only
# roundoff: the identity dominates the weight of E, and at the chain's
# own cutoff the z^2/2 H^2 term, about z^2 of it, is discarded once nt0
# is large (measured on 8 sites at nt0=4000 and cutoff 1e-16: <H> drifted
# by 1.9e-4 over t=1). Every MPS application is truncated at maxdim and
# the chain's cutoff.
function exponential(H,psi,tau,nt0,maxdim,mpomaxdim,cutoff,sites)
	z = tau/nt0 # small step
	iden = identity_mpo(sites) # identity operator
	H2 = apply(H,H;maxdim=mpomaxdim,cutoff=1e-28) # H^2
	E = add(iden,z*H,(0.5*z*z)*H2;alg="directsum")
	E = truncate(E;maxdim=mpomaxdim,cutoff=1e-28)
	for i=1:nt0
		psi = apply_op(E,psi,maxdim,cutoff)
	end
	return psi
end




function overlap(wf1,wf2)
	c = inner(wf1,wf2)
	return c
end



function random_state(sites)
	return random_mps(sites)
end


# The start of every julia_live DMRG-family solve that has no state of its
# own to sweep from (get_gs_dmrg's fresh branch, both generalized solvers,
# the NH-DMRG attempts, each excited state): a random MPS of link
# dimension `linkdims`, the counterpart of mpscpp3's
# default_mps() = randomMPS(sites_,m) and of pyitensor's _default_mps at
# its ramp's first maxdim. random_state above is left as it is, a product
# state, for the consumers that want a random vector rather than a solver
# start (chain.random_state(), the random-witness probes).
#
# They used to share it, and a bond-dimension-1 start is a trap for
# two-site DMRG on a Hamiltonian that couples two sites across a site
# carrying no term (decoupled sublattices; every Thermal_Spin_Chain's
# physical H on the doubled chain): the two-site update never holds a
# coupled pair, sees the far site only through a bond of dimension 1, and
# without noise returns a product state -- -0.5 against an exact -1.0 on
# the Heisenberg chain on the even sites of 6 (2026-09-25 hole hunt 25b,
# finding 7). The non-Hermitian solves take no noise (nhdmrg.jl's note),
# so the start is the only thing that rescues them.
function dmrg_start_state(sites, linkdims)
	return random_mps(sites; linkdims = linkdims)
end


function mpstimesscalar(a,wf1)
	return a*wf1
end


