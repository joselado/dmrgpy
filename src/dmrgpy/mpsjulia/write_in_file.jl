function write_in_file(name::String,a::Float64,mode="w")
  open(name, mode) do io
          write(io, string(a)) # write energy in a file
	  write(io, string("\n")) # write
  end
end


function write_in_file(name::String,a::Complex,mode="w")
  open(name, mode) do io
	  write(io, string(real(a))) # write energy in a file
	  write(io, string("  ")) # write energy in a file
	  write(io, string(imag(a))) # write energy in a file
	  write(io, string("\n")) # write
  end
end


function write_in_file(name::String,a::Int,mode="w")
  open(name, mode) do io
          write(io, string(a)) # write energy in a file
	  write(io, string("\n")) # write
  end
end


function write_in_file(name::String,a::Array,mode="w")
  open(name, mode) do io
	  for i=1:length(a)
		  write(io, string(a[i])) # write energy in a file
		  write(io, "  ") # write energy in a file
	  end
	  write(io, string("\n")) # write
  end
end


