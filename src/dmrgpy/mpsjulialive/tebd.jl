
# TEBD: 2nd-order Trotter (odd/even bond) real-time evolution for strictly
# nearest-neighbor Hamiltonians, run entirely in Julia. Julia counterpart of
# pyitensor/tebd.py (see that module's own docstring for the algorithm's
# full derivation and why it's correct for any parity-conserving fermionic
# Hamiltonian restricted to two adjacent sites) and mpscpp3/tebd.h (the C++
# port). One call per full nt-step trajectory, same design as tdvp.jl's
# quench_tdvp/evolve_and_measure_tdvp.
#
# Unlike tdvp_step (which hands a full Hamiltonian MPO to ITensorMPS's own
# tdvp()), TEBD needs each term's own *local* per-bond operator -- there is
# no way to slice that out of a compiled AutoMPO/MPO any more here than in
# C++ (automaton bookkeeping inflates even a single term's bond dimension;
# see tebd.h's own comment) -- so this file resolves each Hamiltonian term
# into its per-site matrices from scratch (read_terms/term_add!/
# resolve_term below), a direct port of pyitensor/autompo.py's
# HTerm.add()/HTerm.resolve().
#
# Deliberately fed *raw* term lists (mpo.py's text_mpo()/MO2list() shape,
# still using dmrgpy's own "C"/"Cdag" names) rather than
# MultiOperator.to_terms()'s Jordan-Wigner-*predressed* ("A"/"Adag"/"F")
# shape that mpscpp3/pyitensor's own TEBD consume: confirmed directly that
# real ITensors.jl's builtin "Fermion" site type
# (ITensors/src/lib/SiteTypes/src/sitetypes/fermion.jl) only defines
# op("C",..)/op("Cdag",..)/op("F",..) -- no "A"/"Adag" -- unlike dmrgpy's
# own hand-rolled pyitensor/C++ site tables, which do define "A"/"Adag" as
# dmrgpy's own name for the bare (string-free) annihilation/creation
# operator. ITensors.jl's op("C",s) *is* exactly that same bare, string-free
# operator (no JW string baked in -- the string only gets inserted by
# AutoMPO's automaton, or here, by term_add!/resolve_term below), so raw
# "C"/"Cdag" is both the form real ITensors.jl understands natively and the
# one dmrgpy already serializes for the existing MPO path (mpo.py's
# text_mpo/get_MPO) -- reused verbatim, via read_terms() below, rather than
# inventing a second serialization.
#
# is_fermionic_opname/term_add!/resolve_term operate on these raw names, so
# (unlike the pyitensor/C++ ports, where the equivalent carry-tracking logic
# is mostly an inert pass-through over already-JW-dressed input) the carry
# tracking here does the actual Jordan-Wigner threading.

using ITensors
using ITensorMPS
using LinearAlgebra


# Parses mpo.py's text_mpo()/read_operator.jl's read_operator() shared line
# format directly into a term list [(coef,[(name,site),...]),...] instead of
# compiling an AutoMPO -- same parsing loop as read_operator.jl's
# read_operator(), just building a plain list rather than an AutoMPO.
function read_terms(ls::Vector{String})
	nterms = parse(Int, ls[1])
	terms = Vector{Tuple{ComplexF64, Vector{Tuple{String, Int}}}}()
	for i in 1:nterms
		np = parse(Int, ls[2*i])
		l = split(ls[2*i+1], "  ")
		c = parse(Float64, l[1]) + im*parse(Float64, l[2])
		ops = Tuple{String, Int}[]
		for j in 1:np
			opname = String(l[2*j+1])
			ind = parse(Int, l[2*j+2])
			push!(ops, (opname, ind))
		end
		push!(terms, (c, ops))
	end
	return terms
end


# Name-based classification, mirrors ITensor's own isFermionic()
# (mpscpp3/ITensor/itensor/mps/autompo.cc) and pyitensor's
# sites/base.py::is_fermionic(): any operator name starting with 'C' is
# fermionic for Jordan-Wigner purposes.
is_fermionic_opname(name) = !isempty(name) && name[1] == 'C'


# Insertion-sort a term's raw (possibly site-unordered) factor list into
# site-ascending order, one factor at a time, picking up the sign from
# anticommuting a fermionic factor past every fermionic factor already fixed
# at a larger site -- direct port of pyitensor/autompo.py's HTerm.add().
# Mutates `ops` in place; `coefref` is a Ref so the running sign flip
# survives across repeated calls for the same term.
function term_add!(ops, coefref, opname, site)
	idx = 1
	while idx <= length(ops) && ops[idx][2] <= site
		idx += 1
	end
	if idx <= length(ops) && is_fermionic_opname(opname)
		nflip = count(o -> is_fermionic_opname(o[1]), ops[idx:end])
		if isodd(nflip)
			coefref[] *= -1
		end
	end
	insert!(ops, idx, (opname, site))
end


# The (dim,dim) matrix of single-site operator `opname` at site `s`, in the
# convention that plain matrix-vector multiplication reproduces op(opname,s)
# acting on a bare physical-index vector -- same convention/primitive as
# metts.jl's metts_site_operator_matrix.
function site_local_matrix(s, opname)
	O = op(opname, s)
	return Array(O, s', s)
end


# The full per-site list of (dim_i,dim_i) matrices (one per site, 1..n) that
# reconstructs this term's action once Kronecker-multiplied -- direct port
# of pyitensor/autompo.py's HTerm.resolve(): sites the term doesn't touch
# get plain Id, except while an odd number of fermionic operators have been
# "seen" scanning left to right, when they get F instead (threading the
# Jordan-Wigner string); a touched site composes its own factors (in their
# original relative order) and may pick up an extra trailing F of its own
# depending on the parity carried in from its left versus its own.
function resolve_term(sites, ops)
	n = length(sites)
	by_site = Dict{Int, Vector{String}}()
	for (name, site) in ops
		push!(get!(by_site, site, String[]), name)
	end
	out = Vector{Matrix{ComplexF64}}(undef, n)
	carry = false
	for site in 1:n
		s = sites[site]
		names = get(by_site, site, nothing)
		if names === nothing
			out[site] = site_local_matrix(s, carry ? "F" : "Id")
			continue
		end
		is_site_fermionic = isodd(count(is_fermionic_opname, names))
		need_F = (carry != is_site_fermionic)
		mats = [site_local_matrix(s, nm) for nm in names]
		if need_F
			push!(mats, site_local_matrix(s, "F"))
		end
		combined = mats[1]
		for m in mats[2:end]
			combined = combined * m
		end
		out[site] = combined
		carry = (carry != is_site_fermionic)
	end
	return out
end


# The 1-based (lo,hi) site range where `mats` differs from the local
# identity -- i.e. the sites this term actually, non-trivially acts on.
# `nothing` if the term is trivial everywhere (a bare scalar term, e.g. a
# ground-state-energy shift folded in as coef*Identity). Mirrors
# pyitensor/tebd.py's _true_span() -- see its own docstring for why this
# must be checked on the *resolved* matrices, not the raw op-name list.
function true_span(mats)
	touched = Int[]
	for (i, m) in enumerate(mats)
		d = size(m, 1)
		if !isapprox(m, Matrix{ComplexF64}(I, d, d); atol = 1e-8)
			push!(touched, i)
		end
	end
	isempty(touched) && return nothing
	return touched[1], touched[end]
end


# {b => ITensor} for b=1..n-1 (the bond between sites b and b+1), each an
# operator over (sites[b]',sites[b],sites[b+1]',sites[b+1]) -- built
# entirely from ITensor outer products (op(name,s), which already carries
# the correct (s',s) index orientation, and plain ITensor `*` for tensor-
# product placement onto a bond) rather than dense-matrix kron+reshape, so
# there is no manual index-ordering/array-layout bookkeeping to get wrong.
# One-site (onsite) terms are split half-and-half onto each of a site's
# neighboring bonds (full weight at a chain boundary) -- same convention as
# pyitensor/tebd.py's bond_hamiltonians(); a term spanning 3+ sites raises,
# since TEBD only supports strictly nearest-neighbor Hamiltonians.
function bond_hamiltonians(sites, terms)
	n = length(sites)
	if n < 2
		error("bond_hamiltonians: TEBD needs at least 2 sites")
	end
	H = Dict{Int, ITensor}()
	for b in 1:(n - 1)
		H[b] = ITensor(sites[b]', sites[b], sites[b+1]', sites[b+1])
	end
	for (coef, raw_ops) in terms
		ops = Tuple{String, Int}[]
		coefref = Ref(ComplexF64(coef))
		for (name, site) in raw_ops
			term_add!(ops, coefref, name, site)
		end
		c = coefref[]
		mats = resolve_term(sites, ops)
		span = true_span(mats)
		if span === nothing
			H[1] += c * (op("Id", sites[1]) * op("Id", sites[2]))
			continue
		end
		lo, hi = span
		if hi - lo > 1
			error("TEBD requires a nearest-neighbor Hamiltonian; found a " *
				  "term with non-trivial support spanning sites $lo..$hi")
		end
		if hi == lo
			has_left = lo > 1
			has_right = lo < n
			weight = (has_left && has_right) ? 0.5 : 1.0
			local_tensor = ITensor((weight * c) .* mats[lo], sites[lo]', sites[lo])
			if has_left
				b = lo - 1
				H[b] += op("Id", sites[b]) * local_tensor
			end
			if has_right
				b = lo
				H[b] += local_tensor * op("Id", sites[b+1])
			end
		else
			tA = ITensor(mats[lo], sites[lo]', sites[lo])
			tB = ITensor(mats[hi], sites[hi]', sites[hi])
			H[lo] += c * (tA * tB)
		end
	end
	return H
end


# The 2-site gate exp(-i*tau*Hbond) for one bond, via ITensors.jl's own
# tensor exponential (exp(A::ITensor), which auto-splits A's primed/
# unprimed indices into left/right) -- the Julia-native equivalent of
# mpscpp3/tebd.h's ITensor BondGate (also an ITensor-native exponential, not
# a manual expm) and pyitensor/tebd.py's scipy.linalg.expm (there, a manual
# dense exponential, needed because pyitensor has no ITensor-level exp()
# primitive of its own).
function bond_gate(Hbond, tau)
	return exp(-im * tau * Hbond)
end


# Contracts `gate` onto bond (i,i+1) and SVD-truncates, mutating psi in
# place. Canonicalizes psi's orthogonality center to i first (via
# orthogonalize!, the same primitive tdvp_step/entropy.jl's bond_entropy
# use), which is what makes the local SVD truncation-optimal regardless of
# which bond was last touched.
function apply_bond_gate!(psi, i, gate, cutoff, maxdim)
	orthogonalize!(psi, i)
	s_i = siteind(psi, i)
	theta = noprime(gate * psi[i] * psi[i+1], "Site")
	if i > 1
		left_link = commonind(psi[i], psi[i-1])
		U, S, V, _ = svd(theta, left_link, s_i; cutoff = cutoff, maxdim = maxdim)
	else
		U, S, V, _ = svd(theta, s_i; cutoff = cutoff, maxdim = maxdim)
	end
	psi[i] = U
	psi[i+1] = S * V
end


# Precomputes every bond's dt/2 (odd bonds) and dt (even bonds) evolution
# gate once from a static, nearest-neighbor Hamiltonian -- step() only ever
# needs those two, not both gates for every bond (2nd-order Trotter:
# exp(-i dt/2 H_odd) exp(-i dt H_even) exp(-i dt/2 H_odd)).
struct TEBDGates
	odd::Vector{Int}
	even::Vector{Int}
	gates_half::Dict{Int, ITensor}
	gates_full::Dict{Int, ITensor}
end

function build_tebd_gates(sites, terms, dt)
	n = length(sites)
	h_bonds = bond_hamiltonians(sites, terms)
	odd = collect(1:2:(n - 1))
	even = collect(2:2:(n - 1))
	gates_half = Dict(b => bond_gate(h_bonds[b], dt / 2.0) for b in odd)
	gates_full = Dict(b => bond_gate(h_bonds[b], dt) for b in even)
	return TEBDGates(odd, even, gates_half, gates_full)
end

function tebd_step!(psi, gates::TEBDGates, cutoff, maxdim)
	for b in gates.odd
		apply_bond_gate!(psi, b, gates.gates_half[b], cutoff, maxdim)
	end
	for b in gates.even
		apply_bond_gate!(psi, b, gates.gates_full[b], cutoff, maxdim)
	end
	for b in gates.odd
		apply_bond_gate!(psi, b, gates.gates_half[b], cutoff, maxdim)
	end
	return psi
end


# TEBD counterpart of tdvp.jl's quench_tdvp(): `ls` is a read_terms()-shaped
# text encoding of the (already energy-shifted) Hamiltonian's raw term list
# (mpo.py's text_mpo(), *not* an MPO); A1mpo/A2mpo/wf0 stay MPOs/MPS, same
# as quench_tdvp, since only the Hamiltonian needs per-term local gates --
# the observables are still applied via ordinary MPO application
# (apply_clean, reused unchanged from tdvp.jl).
function quench_tebd(ls, sites, A1mpo, A2mpo, wf0, nt, dt, cutoff, maxdim)
	terms = read_terms(ls)
	gates = build_tebd_gates(sites, terms, dt)
	psi1 = apply_clean(A1mpo, wf0, maxdim, cutoff)
	psi2 = apply_clean(A2mpo, wf0, maxdim, cutoff)
	norm0 = sqrt(abs(inner(psi1, psi1)))
	correlator = ComplexF64[]
	for it = 1:nt
		# Measure before evolving, so correlator[k] is C((k-1)*dt) and
		# lines up with timedependent.py's ts=[0,dt,...,(nt-1)*dt] grid --
		# see evolution_dmrg_DC()'s own comment there. Kept in lockstep with
		# the C++/pyitensor/ED copies of this loop.
		push!(correlator, inner(psi2, psi1))
		psi1 = tebd_step!(psi1, gates, cutoff, maxdim)
		psi1 = psi1 * (norm0 / norm(psi1))
	end
	return correlator, psi1
end


# TEBD counterpart of tdvp.jl's evolve_and_measure_tdvp() -- see that
# function's own comment for why renormalization here targets wf's own
# initial norm rather than 1. `wf` is the caller's own MPS object (unlike
# quench_tebd's psi1/psi2, which are always a fresh apply_clean() result),
# so it's copied before any mutation -- same reasoning as tdvp_step's own
# top-of-function copy().
function evolve_and_measure_tebd(ls, sites, Aop, wf, nt, dt, cutoff, maxdim)
	terms = read_terms(ls)
	gates = build_tebd_gates(sites, terms, dt)
	psi = copy(wf)
	norm0 = sqrt(abs(inner(psi, psi)))
	correlator = ComplexF64[]
	for it = 1:nt
		# Measure before evolving, so correlator[k] is C((k-1)*dt) and
		# lines up with timedependent.py's ts=[0,dt,...,(nt-1)*dt] grid --
		# see evolution_dmrg_DC()'s own comment there. Kept in lockstep with
		# the C++/pyitensor/ED copies of this loop.
		push!(correlator, inner(psi, Aop, psi))
		psi = tebd_step!(psi, gates, cutoff, maxdim)
		psi = psi * (norm0 / norm(psi))
	end
	return correlator, psi
end
