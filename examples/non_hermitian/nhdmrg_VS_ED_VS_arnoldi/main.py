# Non-Hermitian DMRG (NH-DMRG, ITensor v3 backend only) cross-checked
# against exact diagonalization and against the MPS Arnoldi route that
# the other backends use for non-Hermitian Hamiltonians.
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
fc = fermionchain.Fermionic_Chain(n) # create the fermion chain
fc.setup_cpp(version=3) # NH-DMRG only exists on the ITensor v3 backend

h = 0
for i in range(n-1):
    h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i] # hopping
for i in range(n):
    h = h + 1j*(-1)**i*fc.Cdag[i]*fc.C[i] # staggered imaginary potential
for i in range(n-1):
    h = h + (fc.N[i]-0.5)*(fc.N[i+1]-0.5) # interaction
fc.set_hamiltonian(h)

# exact reference: diagonalize the full non-Hermitian matrix,
# sorted by real part
es_ed = fc.get_excited(mode="ED",n=6)
print("ED spectrum (lowest real parts):",es_ed[:3])

# NH-DMRG: (energy, left eigenvector, right eigenvector)
e,psil,psir = fc.nhdmrg()
print("NH-DMRG energy:",e)

# the pre-existing Arnoldi route (still the non-Hermitian fallback for
# the v2/pure-Python backends)
from dmrgpy import mpsalgebra
ea,wfa = mpsalgebra.lowest_energy_non_hermitian_arnoldi(fc,h,n=1)
print("Arnoldi energy:",ea[0])

# NH-DMRG reproduces the smallest real part, and lands on an actual
# eigenvalue (either member of a conjugate pair is acceptable)
assert abs(e.real-es_ed[0].real)<1e-6
assert min(abs(e-x) for x in es_ed)<1e-6
# Arnoldi agrees at its own (much looser) accuracy
assert abs(ea[0].real-es_ed[0].real)<5e-2

# (psil, psir) is a genuine biorthogonal eigenpair of H
r = h*psir - e*psir
print("right-eigenvector residual:",abs(r.dot(r))**0.5)
assert abs(r.dot(r))**0.5<1e-6
l = h.get_dagger()*psil - np.conj(e)*psil
print("left-eigenvector residual:",abs(l.dot(l))**0.5)
assert abs(l.dot(l))**0.5<1e-6
assert abs(psil.dot(psir)-1.0)<1e-8 # <psil|psir> = 1
assert abs(psil.aMb(h,psir)/psil.dot(psir)-e)<1e-8 # Rayleigh quotient

print("All checks passed")
