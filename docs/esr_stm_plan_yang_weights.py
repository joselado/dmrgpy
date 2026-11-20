"""Dimer of Yang et al. SM eq. S7: weights of S14 against the plan's DC
weights (P_m - P_n)|<n|S1x|m>|^2 (<n|S1z|n> - <m|S1z|m>), numerically."""
import numpy as np
sx = np.array([[0, 1], [1, 0]])/2.; sy = np.array([[0, -1j], [1j, 0]])/2.
sz = np.array([[1, 0], [0, -1]])/2.; I2 = np.eye(2)
S1 = [np.kron(o, I2) for o in (sx, sy, sz)]; S2 = [np.kron(I2, o) for o in (sx, sy, sz)]
B, J, D, kT = 25.2, 1.0, 0.0, 25.0     # GHz: 0.9 T at g=2, J, D, 1.2 K
for Bt in (0.2, 1.0, 3.0):
    H = (B+Bt)*S1[2] + B*S2[2] + (J+2*D)*S1[2]@S2[2] + (J-D)*(S1[0]@S2[0]+S1[1]@S2[1])
    e, v = np.linalg.eigh(H); p = np.exp(-(e-e[0])/kT); p /= p.sum()
    X = v.conj().T@S1[0]@v; Z = np.real(np.diag(v.conj().T@S1[2]@v))
    xi = np.arctan2(J-D, Bt)                     # tan(xi) = (J-D)/(gamma hbar Btip)
    print("Btip=%.1f GHz  xi=%.3f  sin^2(xi/2)=%.4f" % (Bt, xi, np.sin(xi/2)**2))
    for n in range(4):
        for m in range(4):
            if e[n] > e[m] and abs(X[n, m]) > 1e-9:
                w = abs(X[n, m])**2
                print("   f=%7.3f  S14 P_m|X|^2=%.4f   plan (P_m-P_n)|X|^2 dZ=%.4f   dZ=%+.4f"
                      % (e[n]-e[m], p[m]*w, (p[m]-p[n])*w*(Z[n]-Z[m]), Z[n]-Z[m]))
