using ITensors
using ITensorMPS
import Random

# The two operations of the CLIK-MPS driver (cliksession.py, through
# mpsjulialive/clik.py) that have no counterpart elsewhere in this folder.
# The rest of what the driver needs is mpsalgebra.jl's overlap and
# mpstimesscalar, mpo_ops.jl's aMb_mpo, dagger_mpo and mult_mpo, and
# tdvp.jl's tdvp_step.


# The two random starts of get_krylov_spectrum, drawn one after the other
# from one generator seeded with `seed`, so that a run is reproducible.
# They are ITensorMPS's random circuit states of link dimension m, real,
# as the Gaussian starts of v3 and "python" are, but not the same draws.
function clik_random_starts(sites, m, seed)
	rng = Random.MersenneTwister(seed)
	psi_p = random_mps(rng, sites; linkdims = m)
	psi_u = random_mps(rng, sites; linkdims = m)
	return psi_p, psi_u
end


# The global subspace expansion before the first steps of each branch, the
# counterpart of pyitensor/gse.py, the reference, and of v3's
# Chain::global_subspace_expand(G,psi,order,cutoff,maxdim): the order-1
# Krylov vectors G psi, G^2 psi, ..., each applied at maxdim and at the
# chain cutoff and normalized, enter ITensorMPS's own
# expand(...; alg="orthogonalize"), and a bond that comes out above maxdim
# is cut back to maxdim with no cutoff, as v3 does. A cutoff there would
# drop exactly the directions the expansion added, since they carry zero
# weight in psi (v3's note in chain_session.h).
#
# The cutoff of the Krylov vectors follows "python", where gse.py applies
# them with applyMPO(H, cur, cutoff=cutoff, maxdim=maxdim), at the chain
# cutoff that clik.py hands over. v3 differs: with maxdim > 0 its addBasis
# gives applyMPO only MaxDim (TDVP/basisextension.h), so each vector is
# cut at ITensor's default cutoff of 1e-13 for the DensityMatrix method.
#
# Not tdvp.jl's gse_expand, which serves tevol_method="TDVP_GSE": it hands
# its order to expand(...; alg="global_krylov") as krylovdim, the number
# of Krylov vectors, one more than v3 and "python" take at the same order,
# and it caps each vector at maxlinkdim(psi)+1 instead of maxdim.
# projection_cutoff is the 1e-12 of v3's addBasis and of pyitensor/gse.py,
# below which a bond gets no new direction; ITensorMPS's default is
# sqrt(eps).
function clik_expand(G, psi, order, cutoff, maxdim)
	psi = noprime(copy(psi), "Link")
	orthogonalize!(psi, 1)
	refs = MPS[]
	cur = psi
	for _ in 1:(order - 1)
		cur = apply(G, cur; maxdim = maxdim, cutoff = cutoff,
		            alg = "densitymatrix")
		nrm = norm(cur)
		nrm > 0 || break # G annihilates psi: nothing to add
		cur = (1 / nrm) * cur
		push!(refs, cur)
	end
	isempty(refs) && return psi
	psi = expand(psi, refs; alg = "orthogonalize", cutoff = cutoff,
	             projection_cutoff = 1e-12)
	if maxlinkdim(psi) > maxdim
		truncate!(psi; maxdim = maxdim) # no cutoff, see above
	end
	return psi
end
