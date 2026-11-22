"""Which reference makes the best Gaussian basis of an Anderson impurity,
with U written as U n_up n_dn: the bare level, the level shifted by U/2,
the restricted Hartree level, or the natural orbitals of a solved state.
Same v3 backend, 12 sweeps, cold start for every solve.

python3 reference.py <uniform|wilson> <n> <U> <ed> <maxms> [nsweeps]"""
import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'src'))
import numpy as np, time, json
from dmrgpy import fermionchain, wilsonchain

def uniform(n, V=1.0, t=1.0):
    h = np.zeros((n,n)); h[0,1]=h[1,0]=-V
    for i in range(1,n-1): h[i,i+1]=h[i+1,i]=-t
    return h
def wilson(n, tD=0.2, Lambda=2.0):
    eps, t = wilsonchain.wilson_chain(n-1, Lambda=Lambda)
    h = np.diag(np.concatenate([[0.0], eps])); h[0,1]=h[1,0]=tD
    for k in range(n-2): h[k+1,k+2]=h[k+2,k+1]=t[k]
    return h
def dims(fc, wf):
    s = fc._session
    mpo = max(np.asarray(a).shape[-1] for a in s.mpo_site_arrays(s.hamiltonian_mpo()))
    mps = max(np.asarray(a).shape[-1] for a in s.mps_site_arrays(wf.cpp_handle))
    return int(mpo), int(mps)

case, n, U, ed = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
maxms = [int(x) for x in sys.argv[5].split(",")]
nsweeps = int(sys.argv[6]) if len(sys.argv) > 6 else 12
h = uniform(n) if case == "uniform" else wilson(n)  # h[0,0] = 0

def hamiltonian(fc):
    H = sum(h[i,j]*(fc.Cdagup[i]*fc.Cup[j]+fc.Cdagdn[i]*fc.Cdn[j]) for i in range(n) for j in range(n) if h[i,j] != 0)
    return H + ed*(fc.Nup[0]+fc.Ndn[0]) + U*fc.Nup[0]*fc.Ndn[0]

def filled(level):  # U=0 state at the level, filled below the Fermi level 0
    hh = h.copy(); hh[0,0] = level
    e, u = np.linalg.eigh(hh)
    tol = 1e-12*np.max(np.abs(e))
    f = np.where(e < -tol, 1.0, np.where(e > tol, 0.0, 0.5))
    return (u.conj()*f) @ u.T

def hartree():  # restricted: level ed + U nbar, nbar its own occupation
    lo, hi = 0.0, 1.0
    for _ in range(60):
        nb = 0.5*(lo+hi)
        if nb - filled(ed+U*nb)[0,0].real > 0: hi = nb
        else: lo = nb
    return ed + U*0.5*(lo+hi)

def chain(m):
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    fc.maxm = m; fc.nsweeps = nsweeps
    return fc

out = dict(case=case, n=n, U=U, ed=ed, nsweeps=nsweeps, runs=[])
fname = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'ref_%s_n%d_U%g_ed%g%s.json' % (case, n, U, ed,
                     "" if nsweeps == 12 else "_ns%d" % nsweeps))
def run(name, b, m):
    fc = chain(m)
    H = hamiltonian(fc)
    t0 = time.time()
    H = H if b is None else b.rotate(H)
    fc.set_hamiltonian(H)
    e = fc.gs_energy(); t = time.time()-t0
    wf = fc.get_gs()
    S = [float(wf.get_bond_entropy(i)) for i in range(fc.ns-1)]
    mpo, mps = dims(fc, wf)
    r = dict(basis=name, maxm=m, E=float(e), S=S, mpo=mpo, t=t)
    out["runs"].append(r)
    print("%-14s maxm=%3d E=%.10f Smax=%.3f mpo=%d t=%.1fs" % (name, m, e, max(S), mpo, t), flush=True)
    json.dump(out, open(fname, 'w'))
    return fc

levels = dict(bare=ed, half=ed+U/2, hartree=hartree())
out["levels"] = levels
print("levels", levels, flush=True)
bases = {"site": None}
for name, level in levels.items():
    bases[name] = chain(8).get_gaussian_basis(cm=filled(level), fixed=[0])
for name, b in bases.items():
    for m in maxms:
        fc = run(name, b, m)
# natural orbitals of the hartree-basis state at maxm=32, at three cutoffs
src = run("hartree", bases["hartree"], 32)
t0 = time.time()
bases["hartree"].chain = src
cm = bases["hartree"].original_correlation_matrix(wf=src.get_gs())
cm = 0.5*(cm+cm.conj().T)
out["spin_mixing"] = float(np.max(np.abs(cm[0::2, 1::2])))
out["occupations"] = [float(x) for x in np.linalg.eigvalsh(cm[0::2, 0::2])]
cm[0::2, 1::2] = 0.0; cm[1::2, 0::2] = 0.0  # as natural_basis does
print("spin mixing dropped %.2e, correlation matrix in %.1fs" % (out["spin_mixing"], time.time()-t0), flush=True)
for cut in (1e-8, 1e-5, 1e-3):
    b = src.get_gaussian_basis(cm=cm, fixed=[0], eigval_cutoff=cut)
    print("natural cutoff %g: basis in %.1fs" % (cut, time.time()-t0), flush=True)
    for m in maxms:
        fc = run("natural%g" % cut, b, m)
# a second iteration, from the natural basis (1e-8) at maxm=32
b1 = src.get_gaussian_basis(cm=cm, fixed=[0])
fc = run("natural1e-08", b1, 32)
b1.chain = fc
b2 = b1.natural_basis(wf=fc.get_gs())
for m in maxms:
    run("natural2", b2, m)
run("natural1e-08", b1, 64)
