

n = 100 # maximum number of operators

f = open("ampotk.h","w")

for ni in range(1,n): # number of operators
    if ni==1: f.write("if (numprod==1) {\n")
    else: f.write("else if (numprod=="+str(ni)+") {\n")
    f.write("  string ") # define the strings
    for i in range(ni): # loop over number
        f.write("op"+str(i)) # define
        if i<(ni-1): f.write(",")
    f.write("; \n")
    f.write("  int ") # define the strings
    for i in range(ni): # loop over number
        f.write("i"+str(i)) # define
        if i<(ni-1): f.write(",")
    f.write("; \n")
    # read the data
    f.write("  hfile >> cr >> ci")
    for i in range(ni):
        f.write(" >> op"+str(i)+" >> i"+str(i))
    f.write("; \n")
    f.write("  auto cz = cr + ci*1i;\n")
    f.write("  ampo += cz")
    for i in range(ni): f.write(",op"+str(i)+",i"+str(i))
    f.write(";\n}\n\n")
f.write("else exit(EXIT_FAILURE) ;\n")
f.close()

  
