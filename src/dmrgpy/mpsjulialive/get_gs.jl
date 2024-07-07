



function get_gs_dmrg(H,psi0; nsweeps = 10,cutoff=1e-8,maxm=80 )
  sweeps = Sweeps(nsweeps)
  maxdim!(sweeps, maxm)
  cutoff!(sweeps, cutoff)
  open("/dev/null", "w") do devnull
        redirect_stdout(devnull) do
            return dmrg(H,psi0, sweeps)
        end
    end
#  energy, psi = dmrg(H,psi0, sweeps)
#  return energy,psi
end





