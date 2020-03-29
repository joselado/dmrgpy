function write_in_file(name::String,a::Float64)
  open(name, "w") do io
          write(io, string(a)) # write energy in a file
  end
end


function write_in_file(name::String,a::Complex)
  open(name, "w") do io
	  write(io, string(real(a))) # write energy in a file
	  write(io, string("  ")) # write energy in a file
	  write(io, string(imag(a))) # write energy in a file
  end
end

