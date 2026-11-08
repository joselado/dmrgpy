# Site types and operator names that ITensors.jl does not provide, for
# the dmrgpy site-type codes of sites.jl's get_sites:
#
#   S=3/2, S=2, S=5/2 (codes 4, 5, 6): one generic spin-S builder, basis
#     m = S, S-1, ..., -S, the same order as pyitensor/sites/spin.py and
#     mpscpp3's spin headers.
#   Z2, Z3, Z4 (codes -1, -2, -3): the clock/shift pair of
#     pyitensor/sites/parafermion.py, Sig|k> = |k-1 mod n> and
#     Tau = diag(w^k), w = exp(2 pi i/n). The Z_n strings of
#     parafermionchain.py are explicit Tau factors in the terms, so no
#     ITensors fermion machinery is involved.
#   Names missing on the stock types: N0..N98 projectors on Boson sites,
#     Sp/Sm/ISy and friends on S=1/2 and S=1, Sy/S2/FermiPhase on
#     Electron, and an identity "F" on every non-fermionic type (ITensorMPS
#     puts Op("F",n) on every site a fermionic term crosses, and a spin
#     site between two fermion modes of Mixed_Spin_Fermion_Chain must
#     treat it as the identity, as ITensor v3's SiteSet::op() does).
#
# Every matrix here is in the ITensors.jl convention M[out, in], so it is
# the transpose of the corresponding pyitensor/sites matrix.
# A name is only defined when op(name, index) does not already work, so
# no method ITensors.jl ships is overwritten. Loaded before sites.jl.

using ITensors
using ITensorMPS
using LinearAlgebra

# generic spin-S matrices, twoS = 2S
function dmrgpy_spin_matrices(twoS::Int)
  d = twoS + 1
  S = twoS / 2
  m = [S - k for k in 0:(d - 1)]
  Sz = diagm(0 => m)
  Sp = zeros(Float64, d, d)
  for k in 2:d # <m_k + 1| S+ |m_k>, the state above is k-1
    Sp[k - 1, k] = sqrt(S * (S + 1) - m[k] * (m[k] + 1))
  end
  Sm = Matrix(transpose(Sp))
  Sx = (Sp + Sm) / 2
  Sy = (Sp - Sm) / (2im)
  ISy = (Sp - Sm) / 2 # i*Sy, real
  Id = Matrix{Float64}(I, d, d)
  ops = Dict{String,Matrix}(
    "Sz" => Sz, "Sx" => Sx, "Sy" => Sy, "ISy" => ISy,
    "S+" => Sp, "S-" => Sm, "Sp" => Sp, "Sm" => Sm,
    "Splus" => Sp, "Sminus" => Sm,
    "Sz2" => Sz * Sz, "Sx2" => real(Sx * Sx), "Sy2" => real(Sy * Sy),
    "S2" => S * (S + 1) * Id, "Id" => Id, "F" => Id,
  )
  for k in 1:d
    v = zeros(Float64, d, d); v[k, k] = 1.0
    ops["proj$(k)"] = v
  end
  return ops
end

# Z_n clock/shift matrices
function dmrgpy_zn_matrices(n::Int)
  w = exp(2im * pi / n)
  Sig = zeros(Float64, n, n)
  for k in 1:n
    Sig[mod(k - 2, n) + 1, k] = 1.0 # Sig|k> = |k-1 mod n>
  end
  Tau = diagm(0 => [w^k for k in 0:(n - 1)])
  if n == 2
    Tau = real(Tau)
  end
  Id = Matrix{Float64}(I, n, n)
  ops = Dict{String,Matrix}(
    "Sig" => Sig, "SigDag" => Matrix(transpose(Sig)),
    "Tau" => Tau, "TauDag" => Matrix(adjoint(Tau)),
    "N" => diagm(0 => Float64.(0:(n - 1))), "Id" => Id, "F" => Id,
  )
  for k in 0:(n - 1)
    v = zeros(Float64, n, n); v[k + 1, k + 1] = 1.0
    ops["Proj$(k)"] = v
  end
  return ops
end

# does op(name, s) already work for this site type?
function dmrgpy_op_exists(name::String, tag::String, d::Int)
  s = Index(d, "Site,$tag,n=1")
  try
    op(name, s)
    return true
  catch
    return false
  end
end

# The (tag, name) => matrix table of every operator defined here. The
# method bodies read it, and the methods are defined at top level with
# @eval on the type objects themselves (SiteType{Symbol(tag)} written
# as an expression did not register a method).
const DMRGPY_SITE_OPS = Dict{Tuple{String,String},Matrix}()

for (tag, twoS) in (("S=3/2", 3), ("S=2", 4), ("S=5/2", 5),
                    ("S=1/2", 1), ("S=1", 2))
  for (name, M) in dmrgpy_spin_matrices(twoS)
    DMRGPY_SITE_OPS[(tag, name)] = M
  end
end
for (tag, n) in (("Z2", 2), ("Z3", 3), ("Z4", 4))
  for (name, M) in dmrgpy_zn_matrices(n)
    DMRGPY_SITE_OPS[(tag, name)] = M
  end
end
# Electron, basis Emp, Up, Dn, UpDn
DMRGPY_SITE_OPS[("Electron", "Sy")] = ComplexF64[0 0 0 0; 0 0 -0.5im 0;
                                                 0 0.5im 0 0; 0 0 0 0]
DMRGPY_SITE_OPS[("Electron", "S2")] = diagm(0 => [0.0, 0.75, 0.75, 0.0])
DMRGPY_SITE_OPS[("Electron", "FermiPhase")] = diagm(0 => [1.0, -1.0, -1.0, 1.0])

# the dimension of each site type, and the space of the new types (the
# higher spins with an optional Sz quantum number in 2*Sz units)
const DMRGPY_SITE_DIMS = Dict("S=1/2" => 2, "S=1" => 3, "S=3/2" => 4,
    "S=2" => 5, "S=5/2" => 6, "Z2" => 2, "Z3" => 3, "Z4" => 4,
    "Electron" => 4)
for (tag, twoS) in (("S=3/2", 3), ("S=2", 4), ("S=5/2", 5))
  st = typeof(SiteType(tag))
  @eval function ITensors.space(::$st; conserve_qns=false,
      conserve_sz=conserve_qns, qnname_sz="Sz")
    if conserve_sz
      return [QN(qnname_sz, $twoS - 2k) => 1 for k in 0:$twoS]
    end
    return $(twoS + 1)
  end
end
for (tag, n) in (("Z2", 2), ("Z3", 3), ("Z4", 4))
  st = typeof(SiteType(tag))
  @eval ITensors.space(::$st) = $n
end

# define every operator of the table that op() does not already know
for ((tag, name), M) in DMRGPY_SITE_OPS
  dmrgpy_op_exists(name, tag, DMRGPY_SITE_DIMS[tag]) && continue
  on = typeof(OpName(name)) # the type objects, as OpName"..." gives them
  st = typeof(SiteType(tag))
  @eval ITensors.op(::$on, ::$st) = DMRGPY_SITE_OPS[($tag, $name)]
end

# occupation projectors N0..N98 on Boson sites; Boson is an alias of
# Qudit, whose operators take the local dimension as an argument
for k in 0:98
  on = typeof(OpName("N$k"))
  @eval function ITensors.op(::$on, ::SiteType"Qudit", d::Int)
    M = zeros(Float64, d, d)
    $k < d && (M[$k + 1, $k + 1] = 1.0)
    return M
  end
end
