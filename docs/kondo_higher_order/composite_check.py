"""Usage: python3 composite_check.py 0.61,2.17,4.33 2e-5 0.25 extra
(positive bath levels in meV, the J step, U, and "extra" to include the
U-at-the-tunnel-vertex term that conductance.py drops; see
docs/kondo_higher_order_plan.md).

Check: cumulative spectral weight of X_l^dag = sum_l' psi0^dag_l' tau(l',l)
(tau = sigma.S/2 + U) in the exact ground state of impurity + discrete
p-h symmetric bath with exchange J_here psi0^dag (sigma/2) psi0 . S,
expanded to first order in J, against conductance.py's 2nd+3rd order
structure with F replaced by the discrete bath function."""
import sys, itertools
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
import warnings; warnings.filterwarnings("ignore")
import numpy as np, scipy.sparse as sp
from types import SimpleNamespace
from dmrgpy.kondospectrumtk import conductance as cd

MUB = 5.7883818060e-5
SIG = [np.array([[0,1],[1,0]],complex), np.array([[0,-1j],[1j,0]]), np.diag([1.,-1.]).astype(complex)]
def spin(s):
    m = np.arange(s,-s-1,-1); d=len(m); Sp=np.zeros((d,d),complex)
    for a in range(1,d): Sp[a-1,a]=np.sqrt(s*(s+1)-m[a]*(m[a]+1))
    return [(Sp+Sp.T)/2,(Sp-Sp.T)/(2j),np.diag(m).astype(complex)]

eps_pos = 1e-3*np.array([float(x) for x in sys.argv[1].split(",")])
lev = np.concatenate([-eps_pos[::-1], eps_pos]); w = np.concatenate([np.diff(np.r_[0,eps_pos])[::-1], np.diff(np.r_[0,eps_pos])])
Nb = len(lev); NM = 2*Nb
S = spin(1.); dimp = 3
Himp = -1e-3*S[2]@S[2] + 2*MUB*(1.0*S[0] + 2.0*S[2])
U = float(sys.argv[3]); EXTRA = len(sys.argv) > 4 and sys.argv[4] == "extra"

# Jordan-Wigner fermions, mode index 2*q+s
a = sp.csr_matrix(np.array([[0.,1.],[0.,0.]])); z = sp.csr_matrix(np.diag([1.,-1.])); e2 = sp.identity(2,format='csr')
def cop(j):
    out = sp.identity(1,format='csr')
    for n in range(NM): out = sp.kron(out, z if n<j else (a if n==j else e2), format='csr')
    return out
C = [cop(j) for j in range(NM)]; CD = [c.T.tocsr() for c in C]
If = sp.identity(2**NM, format='csr'); Ii = sp.identity(dimp, format='csr')
kr = lambda F, Imp: sp.kron(F, Imp, format='csr')
Nop = sum(CD[j]@C[j] for j in range(NM))
H0 = kr(sum(lev[q]*(CD[2*q+s]@C[2*q+s]) for q in range(Nb) for s in range(2)), Ii) + kr(If, sp.csr_matrix(Himp))
psi = [sum(np.sqrt(w[q])*C[2*q+s] for q in range(Nb)) for s in range(2)]
psid = [p.T.conj().tocsr() for p in psi]
Vex = sum(kr(psid[s1]@psi[s2], sp.csr_matrix(0.5*sum(SIG[al][s1,s2]*S[al] for al in range(3))))
          for s1 in range(2) for s2 in range(2))
Xd = [sum(kr(psid[lp], sp.csr_matrix(0.5*sum(SIG[al][lp,l]*S[al] for al in range(3)) + (U*np.eye(dimp) if lp==l else 0)))
          for lp in range(2)) for l in range(2)]
nfull = np.array(kr(Nop, Ii).diagonal()).real
N0 = Nb  # half filling
sec0 = np.nonzero(np.abs(nfull-N0)<.5)[0]; sec1 = np.nonzero(np.abs(nfull-N0-1)<.5)[0]

def weights(J, Vs):
    H = (H0 + J*Vex).tocsr()
    H0s = H[sec0][:,sec0].toarray(); H1s = H[sec1][:,sec1].toarray()
    E0s, v0 = np.linalg.eigh(H0s); E1s, v1 = np.linalg.eigh(H1s)
    g = v0[:,0]; E0 = E0s[0]
    tot = np.zeros(len(E1s))
    for l in range(2):
        amp = v1.conj().T @ (Xd[l][sec1][:,sec0] @ g)
        tot += np.abs(amp)**2
    return np.array([2*np.pi*2*np.sum(tot[(E1s-E0)<V]) for V in Vs]), E1s-E0

# predicted
e, v = np.linalg.eigh(Himp); e = e-e[0]
Se = [v.conj().T@x@v for x in S]
ks = SimpleNamespace(e=e, p=np.r_[1.,0,0], dim=3, T=0., kB=8.617e-5, Sx=Se[0], Sy=Se[1], Sz=Se[2])
M2 = cd._spin_matrix_elements_squared(ks)
coeff = np.imag(cd._triple_product_coefficients(ks, np.array([0])))[0]/2.  # [f,m]
Xi = np.stack(Se, -1)
loop = np.real(np.einsum('mk,mk->m', Xi[0], Xi[:,0]))
def Fdisc(x, excl=None):
    t = w[Nb//2:]/(eps_pos - x)
    if excl is not None: t = np.where(np.abs(eps_pos-x)<1e-12, 0., t)
    return t.sum()
def predicted(Vs):
    out0 = np.zeros(len(Vs)); out1 = np.zeros(len(Vs))
    for qi, (eq, wq) in enumerate(zip(eps_pos, w[Nb//2:])):
        for f in range(3):
            ev = eq + e[f]  # on-shell bias for this final state
            sel = ev < np.asarray(Vs)
            z0 = 2*np.pi*(M2[f,0] + (4*U**2 if f==0 else 0.))
            k3 = 0.
            for m in range(3):
                xd = ev - e[m]; xx = ev - (e[f]-e[m])
                Fd = Fdisc(xd, excl=True); Fx = Fdisc(-xx)
                k3 += coeff[f,m]*(Fd + Fx)
                if f == 0: k3 += U*loop[m]*(Fd - Fx)
            if EXTRA: k3 += U*M2[f,0]*(Fdisc(ev, excl=True) - Fdisc(-ev))  # U at the tunnel vertex, spin at the exchange vertex, any f
            out0 += sel*wq*z0; out1 += sel*wq*4*np.pi*k3
    return out0, out1

_, ex = weights(0., [0.])
thr = np.unique(np.round(ex[ex<12e-3], 9))
Vs = 0.5*(thr[:-1]+thr[1:])
Vs = Vs[np.diff(thr) > 2e-5]
dJ = float(sys.argv[2])
Wp, _ = weights(+dJ, Vs); Wm, _ = weights(-dJ, Vs); W0, _ = weights(0., Vs)
p0, p1 = predicted(Vs)
dW = (Wp-Wm)/(2*dJ)   # d/dJ_here ; code Jrho_s = -J_here
print("zeroth order max diff:", np.max(np.abs(W0-p0)), " scale", np.max(np.abs(p0)))
print("levels",sys.argv[1],"dJ",dJ,"U",U,"extra term",EXTRA); print("first order: exact dW/dJ  vs  -pred"); print("max |diff|", np.max(np.abs(dW+p1)), " scale", np.max(np.abs(p1)))
for V, a1, b1 in zip(Vs, dW, -p1): print(f"V={V*1e3:7.3f} meV  exact {a1: .6f}  pred {b1: .6f}  diff {a1-b1: .2e}")
