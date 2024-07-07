
function get_sites(ls::Vector{String})
  n = parse(Int,ls[1]) # number of sites
  sites = []
  for i=1:n
    ni = parse(Int,ls[i+1]) # index of the site
    if ni==0 
	site = Index(2,"Site,Fermion,n=$i")
    end 
    if ni==2 
	site = Index(2,"Site,S=1/2,n=$i")
    end 
    if ni==3 
	site = Index(3,"Site,S=1,n=$i")
    end
    push!(sites,site)
  end
#  sites = [s for s in sites]
  return sites
end

