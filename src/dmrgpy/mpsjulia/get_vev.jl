
function get_vev()
    sites = get_sites()
    wf = get_gs(sites) # get the ground state wavefunction
    Am = read_operator("vev_multioperator.in") # read operator
    A = MPO(Am,sites)
    out = inner(wf,A,wf)
    write_in_file("VEV.OUT",out)
#    print(out)
end

