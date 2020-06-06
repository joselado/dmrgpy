

function save_mps(filename,psi)
	serialize(filename,psi)
end


function load_mps(filename)
  return deserialize(filename)
end

