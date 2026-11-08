# Conserved-sector mode on julia_live: a QN-carrying site set built from
# the same site codes as sites.jl's dense one, the in-sector random start,
# and the exact promotion of a QN state back to the dense site set.
#
# Only the quantum numbers the caller asked for are conserved: an Electron
# site with conserve_nf alone leaves Sz free. ITensors uses the same units
# as dmrgpy's set_conserved_sector, Sz in 2*Sz (QN("Sz",+1) for "Up" on
# S=1/2), Nf and Number as counts. The basis states are in the same order
# on the QN and the dense index of every site type below (S=1/2 Up,Dn;
# S=1 Up,Z0,Dn; Fermion Emp,Occ; Electron Emp,Up,Dn,UpDn; Boson 0..d-1),
# which is what makes the promotion a relabelling and nothing else.

function get_sites_qn(ls::Vector{String}, conserve_nf::Bool,
        conserve_sz::Bool, conserve_nb::Bool)
  n = parse(Int,ls[1]) # number of sites
  sites = Any[]
  for i=1:n
    ni = parse(Int,ls[i+1]) # site code
    if ni==0
      site = siteind("Fermion",i; conserve_nf=conserve_nf,
                     conserve_nfparity=false)
    elseif ni==1
      site = siteind("Electron",i; conserve_nf=conserve_nf,
                     conserve_sz=conserve_sz, conserve_nfparity=false)
    elseif ni==2
      site = siteind("S=1/2",i; conserve_sz=conserve_sz)
    elseif ni==3
      site = siteind("S=1",i; conserve_sz=conserve_sz)
    elseif 100 < ni < 200
      site = siteind("Boson",i; dim=ni-100, conserve_number=conserve_nb,
                     qnname_number="Nb")
    else
      error("conserved-sector mode on julia_live has no QN site for site code $ni")
    end
    if !hasqns(site)
      error("site $i (code $ni) conserves none of the requested quantum numbers")
    end
    push!(sites,site)
  end
  # concretely typed, Vector{Index{Vector{Pair{QN,Int}}}}, unlike the
  # dense Vector{Index} of sites.jl: the method of applyoperator below is
  # chosen by that type, so QN site sets reach it and dense ones never do
  return identity.(sites)
end


# MPO times MPS on QN sites. mpsalgebra.jl's applyoperator removes the
# prime of contract()'s output by a second contraction against an identity
# MPO, which QN indices refuse (both primed indices point out). apply()
# does contract() and noprime() in one step. It contracts exactly and then
# truncates ("naive"), since the "densitymatrix" algorithm raises a
# BoundsError on QN blocks that a product of a few sectors leaves empty.
# A result that is zero, exactly (N on the empty state) or to rounding
# (the total Sz on an Sz=0 state), is returned untruncated: truncate()
# raises on it, with no blocks at all or with a norm that is the square
# root of -4e-17. The untruncated product is exact, only larger.
function applyoperator(sites::AbstractVector{<:Index{<:AbstractVector}},
        A, psi0, maxdim, cutoff; alg = "naive")
  psi1 = apply(A, psi0; alg = "naive", truncate = false)
  try
    return truncate(psi1; maxdim = maxdim, cutoff = cutoff)
  catch e
    (e isa DomainError || e isa BoundsError || e isa ArgumentError) || rethrow()
    return psi1
  end
end


# random MPS of link dimension linkdims inside the sector of the product
# state `states` (1-based basis-state numbers, one per site)
function sector_start_state(sites, states::Vector{Int}, linkdims::Int)
  return random_mps(sites, states; linkdims = linkdims)
end


# true when an MPS lives on QN-carrying site indices
function mps_hasqns(psi)
  return any(hasqns, siteinds(psi))
end


# exact conversion of a QN MPS to the dense site set `dense_sites`
function promote_dense_mps(psi, dense_sites)
  phi = dense(psi)
  return replace_siteinds(phi, dense_sites)
end
