using ITensors

function read_operator(name::String)
	ls = readlines(name) # lines in the file
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
#		print(out,"\n") # print
		ampo += out # Add to the MPO
	end
	print(ampo)
	return ampo # return MPO
end


function read_mpo(name::String)
	sites = get_sites()
	ampo = read_operator(name)
	return MPO(ampo,sites)
end

function identity_mpo()
  ampo = AutoMPO()
  ampo +=     ("Id",1)
  sites = get_sites()
  return MPO(ampo,sites)
end


#read_operator("hamiltonian.in")

