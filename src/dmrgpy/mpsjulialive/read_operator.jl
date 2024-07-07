using ITensors

function read_operator(ls::Vector{String})
	nterms = parse(Int,ls[1]) # number of terms in the operator
	ampo = AutoMPO() # create AMPO
	for i=1:nterms # loop over terms
		out = [] # create empy list
		np = parse(Int,ls[2*i]) # number of terms
		l = split(ls[2*i+1],"  ") # operators for this term
		c = parse(Float64,l[1]) + im*parse(Float64,l[2])  # coupling
		push!(out,c) # add this term
		for j=1:np # loop over terms in the product
			op = String(l[2*j+1]) # name of the operator
			ind = parse(Int,l[2*j+2]) # name of the operator
			push!(out,op)
			push!(out,ind)
		end
		out = tuple(out[1:length(out)]...) # convert to tuple
		ampo += out # Add to the MPO
	end
	return ampo # return MPO
end


function toMPO(ls,sites)
	ampo = read_operator(ls)
	return MPO(ampo,sites)
end

function identity_mpo(sites)
  ampo = AutoMPO()
  ampo +=     ("Id",1)
  return MPO(ampo,sites)
end

