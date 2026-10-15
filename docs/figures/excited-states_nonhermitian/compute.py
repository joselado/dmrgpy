# Many-body skin effect of the interacting Hatano-Nelson chain, from NH-DMRG.
# With open ends H = S H0 S^-1, with H0 the chain with symmetric hopping
# t' = sqrt((t+g)(t-g)) and S = prod_j r^(j n_j), r = t'/(t+g): this is the
# exact reference at n sites, and ED is the reference at n_ed sites.
import numpy as np
from dmrgpy import fermionchain

t, V = 1.0, 1.0
gammas = [0.0, 0.1, 0.2, 0.3, 0.4]
n, n_ed = 16, 10


def chain(n, tl, tr):  # tl on c^dag_i c_i+1, tr on c^dag_i+1 c_i
    fc = fermionchain.Fermionic_Chain(n)
    fc.set_hamiltonian(sum(tl*fc.Cdag[i]*fc.C[i+1] + tr*fc.Cdag[i+1]*fc.C[i]
                           + V*(fc.N[i]-0.5)*(fc.N[i+1]-0.5)
                           for i in range(n-1)))
    fc.maxm, fc.nsweeps, fc.cutoff = 80, 20, 1e-12
    return fc


def density(fc, wf, **kw):  # right-vector density <n_i>_R
    norm = fc.aMb(wf, fc.Id, wf, **kw).real
    return [fc.aMb(wf, fc.N[i], wf, **kw).real/norm for i in range(len(fc.N))]


def nh_solve(n, g):
    fc = chain(n, t+g, t-g)
    e, psil, psir = fc.nhdmrg(tol=1e-5, ntries=3)
    res = fc.applyoperator(fc.hamiltonian, psir) - e*psir
    resid = np.sqrt(abs(res.overlap(res))/abs(psir.overlap(psir)))/abs(e)
    return e.real, density(fc, psir), resid


def reference(n, g):
    tp = np.sqrt((t+g)*(t-g))
    fh = chain(n, tp, tp)
    e0, w = fh.gs_energy(), fh.get_gs()
    for j in range(n):  # S, centred to keep the numbers of order one
        w = fh.applyoperator(fh.Id + ((tp/(t+g))**(j-(n-1)/2.)-1)*fh.N[j], w)
    return e0, density(fh, w)


dens, dens_ref, dens_small, dens_ed, converged = [], [], [], [], []
for g in gammas:
    dens_small.append(nh_solve(n_ed, g)[1])
    fe = chain(n_ed, t+g, t-g)
    _, ws = fe.get_excited_states(mode="ED", n=1)
    dens_ed.append(density(fe, ws[0], mode="ED"))
    e, d, resid = nh_solve(n, g)
    e0, dref = reference(n, g)
    dens.append(d)
    dens_ref.append(dref)
    converged.append(abs(e - e0) < 1e-6 and resid < 1e-4)  # else not plotted

np.savez("data.npz", gamma=gammas, n=n, n_ed=n_ed, dens=dens,
         dens_ref=dens_ref, dens_small=dens_small, dens_ed=dens_ed,
         converged=converged)
