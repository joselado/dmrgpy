"""The Floquet steady state of the driven Lindblad equation, the reference
of the rate equation on small systems (get_esr_spectrum(method="floquet")).

The Lindbladian is the one the rate equation is the secular limit of: in
the eigenbasis of H, jumps sqrt(Gamma_{n->m}) |m><n| with the Korringa
rates of rates.py and pure dephasing sqrt(1/Tphi) |n><n| on every state, so
that a coherence rho_nm decays at (Gamma_n + Gamma_m)/2 + 1/Tphi, as in the
rate equation. The drive Omega V cos(w t) enters without the rotating-wave
and without the secular approximation. With rho(t) = sum_k rho_k e^{-ikwt},
|k| <= nharm,

    (L0 + i k w) rho_k + L1 (rho_{k-1} + rho_{k+1}) = 0,
    L1 = -i (Omega/2) [V, .],   Tr rho_0 = 1,

solved as one sparse linear system per frequency. The readout is the same
as lineshape.py: Tr D (rho_0 - rho_eq) and Re[e^{i phase} Tr D rho_1].
Its size is (2 nharm + 1) dim^2, so it is kept to dim <= MAX_DIM.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

MAX_DIM = 64


def signal(lines, f, vrf_over_vdc=0., phase=0., nharm=2):
    """(dc, homodyne) from the Floquet solution, for the Lines object of
    lineshape.py (same rates, populations, drive and units)"""
    spec = lines.spec
    d = spec.dim
    if d > MAX_DIM:
        raise ValueError("get_esr_spectrum(method=\"floquet\"): %d states are "
                         "kept, more than the %d this reference is made for; "
                         "use method=\"rates\"" % (d, MAX_DIM))
    nharm = int(nharm)
    if nharm < 1:
        raise ValueError("nharm must be at least 1")
    I = sp.identity(d, format="csr")
    # L0 in row-major vectorization, index a*d + b for rho_ab
    wa = lines.w_n
    diag = (-1j*(wa[:, None] - wa[None, :])
            - 0.5*(lines.out[:, None] + lines.out[None, :])
            - lines.gphi*(1. - np.eye(d)))
    L0 = sp.diags(diag.ravel(), format="lil")
    Gam = lines.Gamma
    for m in range(d):
        for n in range(d):
            if m != n and Gam[m, n] != 0.:
                L0[m*d + m, n*d + n] += Gam[m, n]
    L0 = sp.csr_matrix(L0)
    V = sp.csr_matrix(np.asarray(spec.V, dtype=complex))
    L1 = -0.5j*lines.Omega*(sp.kron(V, I) - sp.kron(I, V.T))
    D = np.asarray(spec.D, dtype=complex)
    nk = 2*nharm + 1
    N = d*d
    # one population row of the k = 0 block per block of the relaxation is
    # replaced by the weight of that block, which the drive does not change
    # (lineshape.Lines has checked that it connects no two blocks)
    constraints = []
    for blk in np.unique(lines.labels):
        idx = np.flatnonzero(lines.labels == blk)
        row = np.zeros(N, dtype=complex)
        row[[a*d + a for a in idx]] = 1.
        constraints.append((idx[0]*d + idx[0], row, np.sum(lines.p0[idx])))
    rho_eq = np.diag(lines.p0)
    dc = np.zeros(len(f))
    hom = np.zeros(len(f))
    ph = np.exp(1j*phase)
    for i, fi in enumerate(np.asarray(f, dtype=float)):
        w = 2.*np.pi*fi
        blocks = [[None]*nk for _ in range(nk)]
        for j in range(nk):
            k = j - nharm
            blocks[j][j] = L0 + 1j*k*w*sp.identity(N, format="csr")
            if j > 0: blocks[j][j-1] = L1
            if j < nk - 1: blocks[j][j+1] = L1
        A = sp.lil_matrix(sp.bmat(blocks, format="csr"))
        r0 = nharm*N # first row of the k = 0 block
        b = np.zeros(nk*N, dtype=complex)
        for row0, row, weight in constraints:
            A[r0 + row0, :] = 0.
            A[r0 + row0, r0:r0 + N] = row
            b[r0 + row0] = weight
        x = sla.spsolve(sp.csr_matrix(A), b)
        rho0 = x[r0:r0 + N].reshape(d, d)
        rho1 = x[r0 + N:r0 + 2*N].reshape(d, d) # coefficient of e^{-iwt}
        dc[i] = np.real(np.trace(D @ (rho0 - rho_eq)))
        if vrf_over_vdc != 0.:
            hom[i] = vrf_over_vdc*np.real(ph*np.trace(D @ rho1))
    return dc, hom
