# Non-Hermitian DMRG (NH-DMRG), cross-checked on every backend that
# implements it (ITensor v3 C++, ITensor v2 C++, pure Python) against
# exact diagonalization and against the MPS Arnoldi route that was
# previously the only non-Hermitian option.
#
# The model is an interacting fermionic chain with a staggered imaginary
# potential (the same one as ../non_hermitian_chain): its spectrum is
# complex, and the "ground state" convention throughout dmrgpy is the
# eigenvalue with the smallest real part. NH-DMRG (a port of
# ITensorNHDMRG.jl's onesided+fidelity algorithm, see
# src/dmrgpy/nhdmrg.py) returns that eigenvalue together with the
# biorthogonal left/right eigenvector pair (<psil|, |psir>), normalized
# so that <psil|psir> = 1.

# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import fermionchain

n = 6

def build_chain(version):
    fc = fermionchain.Fermionic_Chain(n) # create the fermion chain
    if version=="python": fc.setup_python()
    else: fc.setup_cpp(version=version)
    h = 0
    for i in range(n-1):
        h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i] # hopping
    for i in range(n):
        h = h + 1j*(-1)**i*fc.Cdag[i]*fc.C[i] # staggered imaginary potential
    for i in range(n-1):
        h = h + (fc.N[i]-0.5)*(fc.N[i+1]-0.5) # interaction
    fc.set_hamiltonian(h)
    return fc,h

# exact reference: diagonalize the full non-Hermitian matrix,
# sorted by real part
fc,h = build_chain(3)
es_ed = fc.get_excited(mode="ED",n=6)
print("ED spectrum (lowest real parts):",es_ed[:3])

for version in [3,2,"python"]:
    fc,h = build_chain(version)
    # NH-DMRG: (energy, left eigenvector, right eigenvector)
    e,psil,psir = fc.nhdmrg()
    print("NH-DMRG (itensor_version="+str(version)+") energy:",e)
    # NH-DMRG reproduces the smallest real part, and lands on an actual
    # eigenvalue (either member of a conjugate pair is acceptable)
    assert abs(e.real-es_ed[0].real)<1e-6
    assert min(abs(e-x) for x in es_ed)<1e-6
    # (psil, psir) is a genuine biorthogonal eigenpair of H
    r = h*psir - e*psir
    assert abs(r.dot(r))**0.5<1e-3
    l = h.get_dagger()*psil - np.conj(e)*psil
    assert abs(l.dot(l))**0.5<1e-3
    assert abs(psil.dot(psir)-1.0)<1e-8 # <psil|psir> = 1
    assert abs(psil.aMb(h,psir)/psil.dot(psir)-e)<1e-8 # Rayleigh quotient

# the pre-existing Arnoldi route, for comparison (typically several
# orders of magnitude less accurate than NH-DMRG at comparable cost)
from dmrgpy import mpsalgebra
fc,h = build_chain(3)
ea,wfa = mpsalgebra.lowest_energy_non_hermitian_arnoldi(fc,h,n=1)
print("Arnoldi energy:",ea[0])
assert abs(ea[0].real-es_ed[0].real)<5e-2

print("All checks passed")
