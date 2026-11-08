
# One dense site index per dmrgpy site-type code, the convention of
# pyitensor/sites/siteset.py and mpscpp3/get_sites.h. The tag after "Site"
# is the ITensors.jl SiteType name, so op(name, s) resolves on it; the
# types ITensors.jl does not ship (S=3/2, S=2, S=5/2, Z2, Z3, Z4) are
# defined in sites_custom.jl.
function dmrgpy_site_index(code::Int, i::Int)
  code == 0 && return Index(2, "Site,Fermion,n=$i")
  code == 1 && return Index(4, "Site,Electron,n=$i") # Emp, Up, Dn, UpDn
  code == 2 && return Index(2, "Site,S=1/2,n=$i")
  code == 3 && return Index(3, "Site,S=1,n=$i")
  code == 4 && return Index(4, "Site,S=3/2,n=$i")
  code == 5 && return Index(5, "Site,S=2,n=$i")
  code == 6 && return Index(6, "Site,S=5/2,n=$i")
  code == -1 && return Index(2, "Site,Z2,n=$i")
  code == -2 && return Index(3, "Site,Z3,n=$i")
  code == -3 && return Index(4, "Site,Z4,n=$i")
  if 100 < code < 200 # boson, local dimension code-100
    return Index(code - 100, "Site,Boson,n=$i")
  end
  error("julia_live has no site type for dmrgpy site-type code $code " *
        "(site $i); the known codes are 0, 1, 2, 3, 4, 5, 6, -1, -2, -3 " *
        "and 101..199 for bosons")
end

function get_sites(ls::Vector{String})
  n = parse(Int, ls[1]) # number of sites
  sites = Index[]
  for i in 1:n
    push!(sites, dmrgpy_site_index(parse(Int, ls[i + 1]), i))
  end
  return sites
end
