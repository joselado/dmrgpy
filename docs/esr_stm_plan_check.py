"""Check the secular weak-drive ESR formula against the exact O(Omega^2)
Floquet steady state of a driven Lindblad equation, on two spins 1/2.
Units: GHz (energies as frequencies), T in GHz too."""
import numpy as np

sx = np.array([[0, 1], [1, 0]])/2.
sy = np.array([[0, -1j], [1j, 0]])/2.
sz = np.array([[1, 0], [0, -1]])/2.
I2 = np.eye(2)
def site(o, k): return np.kron(o, I2) if k == 0 else np.kron(I2, o)
S = [[site(o, k) for o in (sx, sy, sz)] for k in range(2)]

def model(J, B, Btip, theta):
    H = sum(J*S[0][a] @ S[1][a] for a in range(3))
    H = H + B*(S[0][2] + S[1][2]) + Btip*S[0][2]
    m = np.array([np.sin(theta), 0., np.cos(theta)])
    V = sum(m[a]*S[0][a] for a in range(3))  # drive = tip-site spin along m
    return H, V, V.copy()                     # detection D = drive

def lindbladian(H, Ls):
    d = H.shape[0]; I = np.eye(d)
    L = -1j*(np.kron(H, I) - np.kron(I, H.T))  # row-major vec: rho -> rho.reshape(-1)
    for A in Ls:
        AdA = A.conj().T @ A
        L = L + np.kron(A, A.conj()) - 0.5*np.kron(AdA, I) - 0.5*np.kron(I, AdA.T)
    return L

def rates(e, kT, G0, gphi, vs):
    """Eigenbasis jumps |m><n| with detailed balance, and pure dephasing."""
    d = len(e); Ls = []; R = np.zeros((d, d))
    for n in range(d):
        for m in range(d):
            if n == m: continue
            de = e[n] - e[m]
            g = G0 if de > 0 else G0*np.exp(de/kT)  # n -> m
            R[m, n] = g
            P = np.outer(vs[:, m], vs[:, n].conj())
            Ls.append(np.sqrt(g)*P)
    for n in range(d):
        P = np.outer(vs[:, n], vs[:, n].conj())
        Ls.append(np.sqrt(2*gphi)*P)
    R -= np.diag(R.sum(axis=0))
    return Ls, R

def floquet(H, V, D, Ls, w, Om, rho0):
    d = H.shape[0]
    L0 = lindbladian(H, Ls)
    comm = lambda X: (np.kron(X, np.eye(d)) - np.kron(np.eye(d), X.T))
    src = 1j*(Om/2)*comm(V) @ rho0.reshape(-1)
    rp = np.linalg.solve(L0 + 1j*w*np.eye(d*d), src).reshape(d, d)
    rm = rp.conj().T
    rhs = 1j*(Om/2)*(comm(V) @ (rp + rm).reshape(-1))
    A = np.vstack([L0, np.eye(d).reshape(1, -1)])  # trace constraint
    b = np.concatenate([rhs, [0.]])
    drho = np.linalg.lstsq(A, b, rcond=None)[0].reshape(d, d)
    return np.trace(D @ drho).real, np.trace(D @ rp)

def secular(e, vs, V, D, R, p0, G0, gphi, w, Om):
    Vm = vs.conj().T @ V @ vs; Dm = vs.conj().T @ D @ vs
    d = len(e); W = np.zeros((d, d)); X = 0j
    out = -np.diag(R)  # total out rate of each state
    for n in range(d):
        for m in range(d):
            if e[n] <= e[m] + 1e-12: continue
            wnm = e[n] - e[m]
            g = 0.5*(out[n] + out[m]) + 2*gphi
            lor = g/((w - wnm)**2 + g**2)
            W[n, m] = W[m, n] = 0.5*Om**2*abs(Vm[n, m])**2*lor
            X += (Om/2)*Vm[n, m]*Dm[m, n]*(p0[m] - p0[n])/(w - wnm + 1j*g)
    Wg = W - np.diag(W.sum(axis=0))
    dp = -np.linalg.lstsq(np.vstack([R, np.ones((1, d))]),
                          np.concatenate([Wg @ p0, [0.]]), rcond=None)[0]
    return np.real(np.sum(dp*np.diag(Dm))), X

kT, G0, gphi, Om = 25., 0.004, 0.002, 1e-4   # 1.2 K; T1 ~ 250 ns; weak drive
for Btip in (0.05, 0.5, 2.0):
    H, V, D = model(J=1.0, B=20., Btip=Btip, theta=0.6)
    e, vs = np.linalg.eigh(H)
    p0 = np.exp(-(e - e[0])/kT); p0 /= p0.sum()
    rho0 = vs @ np.diag(p0) @ vs.conj().T
    Ls, R = rates(e, kT, G0, gphi, vs)
    Dm = vs.conj().T @ D @ vs
    print("Btip=%.2f GHz  levels" % Btip, np.round(e - e[0], 4))
    for (n, m) in ((1, 0), (2, 1), (3, 2), (3, 1), (2, 0)):
        wnm = e[n] - e[m]
        dD = (Dm[n, n] - Dm[m, m]).real
        for dw in (-0.003, 0., 0.003):
            w = wnm + dw
            fd, fx = floquet(H, V, D, Ls, w, Om, rho0)
            sd, sx_ = secular(e, vs, V, D, R, p0, G0, gphi, w, Om)
            print("  %d<-%d f=%.4f dD=%+.3f  DC exact %+.4e secular %+.4e | "
                  "hom exact %+.4e%+.4ej secular %+.4e%+.4ej"
                  % (n, m, w, dD, fd/Om**2, sd/Om**2, fx.real/Om, fx.imag/Om,
                     sx_.real/Om, sx_.imag/Om))
