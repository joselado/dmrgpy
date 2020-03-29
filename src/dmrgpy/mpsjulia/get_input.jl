function get_input_string(name::String,default::String="")
  out = ""
  open("tasks.in") do file
      for ln in eachline(file)
	  l = replace(ln," "=>"")
	  l = split(l,"=")
	  if name==l[1]
		  out = l[2]
          end
      end
  end
  out = string(out)
  if out==""
	  return default
  end
  return out
  end 


get_string(name::String) = get_input_string(name)
get_value(name::String) = parse(Int,get_input_string(name))
get_float(name::String) = parse(Float64,get_input_string(name,"0.0"))
get_int(name::String) = parse(Int,get_input_string(name,"0"))
get_bool(name::String) = parse(Bool,get_input_string(name,"false"))


 

#println(get_input_string("maxm"))
#a = get_value("maxm")
#print(a)
