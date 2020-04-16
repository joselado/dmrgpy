
function get_vev()
    sites = get_sites()
    wf = get_gs(sites) # get the ground state wavefunction
    Am = read_operator("vev_multioperator.in") # read operator
    A = MPO(Am,sites)
    out = inner(wf,A,wf)
    write_in_file("VEV.OUT",out)
#    print(out)
end


function get_many_vev()
    sites = get_sites() # get the sites
    wf = get_gs(sites) # get the ground state wavefunction
    nvev = get_int("num_vev") # number of vev
    write_in_file("VEV.OUT","","w")
    for i=1:nvev
	name = string("vev_multioperator_",i,"_.in") # name of the file
        Am = read_operator(name) # read operator
        A = MPO(Am,sites) # get the operator
        out = inner(wf,A,wf)
        write_in_file("VEV.OUT",out,"a")
    end
#    print(out)
end



