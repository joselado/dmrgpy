# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import time
import numpy as np
import matplotlib.pyplot as plt

from dmrgpy import spinchain, fermionchain
from dmrgpy import mpsalgebra
from dmrgpy.multioperator import MultiOperator, MO2matrix

# Hermiticity of an operator can be settled in two ways in dmrgpy, and
# this script sweeps chain length to compare what each one costs.
#
# The symbolic route puts H-H^dagger in canonical form
# (multioperatortk/canonical.py): every term's factors sorted by site,
# with the sign of the fermionic reordering, and equal terms collected.
# For an ordinary Hamiltonian every term cancels against its own
# reversed partner, which is a proof, so nothing numerical runs at all.
# The route is one-sided: True means proven, False only means the
# canonical form did not collapse, which a Hermitian operator can
# survive when its Hermiticity rests on something the operator names do
# not carry, a same-site identity (Sx Sx = 1/4 on a spin-1/2 site) or
# two factors on one site commuting without both being diagonal.
#
# The numerical route is the fallback for exactly that case: apply
# H-H^dagger to a random low-bond-dimension state and look at the norm.
# Many_Body_Chain.is_hermitian() chains the two, which is what
# gs_energy() and exponential() gate on.


def heisenberg(n):
    """Uniform S=1/2 Heisenberg chain and its Hamiltonian"""
    sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
    h = 0
    for i in range(n-1): h = h + sc.SS(i,i+1)
    sc.set_hamiltonian(h)
    return sc,h


def time_symbolic(sc,h,nrep=5):
    """Cost of the canonical-form proof alone"""
    ts = []
    for i in range(nrep):
        t0 = time.time() ; out = h.is_hermitian() ; ts.append(time.time()-t0)
    return min(ts),out


def time_numerical(sc,h,nrep=5):
    """Cost of the random-witness probe alone, with the proof disabled"""
    old = MultiOperator.is_hermitian
    MultiOperator.is_hermitian = lambda self: False # force the fallback
    try:
        ts = []
        for i in range(nrep):
            t0 = time.time() ; out = mpsalgebra.is_hermitian(sc,h)
            ts.append(time.time()-t0)
        return min(ts),out
    finally:
        MultiOperator.is_hermitian = old


ns = [4,8,12,16,24,32,40]
tsym,tnum = [],[]
for n in ns:
    sc,h = heisenberg(n)
    ta,ra = time_symbolic(sc,h)
    tb,rb = time_numerical(sc,h)
    tsym.append(ta) ; tnum.append(tb)
    print("n = %2d   symbolic %.5f s (%s)   numerical %.5f s (%s)"%(n,ta,ra,tb,rb))
    assert ra and rb # both have to see that a Heisenberg chain is Hermitian

# The canonical form is a rewrite, not an approximation: check directly
# that it builds the same matrix, on a fermionic operator where the
# reordering actually carries signs
fc = fermionchain.Fermionic_Chain(4)
hf = 0
for i in range(3): hf = hf + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
fc.set_hamiltonian(hf)
obj = fc.get_ED_obj()
probes = {
        "hopping written backwards": fc.Cdag[2]*fc.C[0],
        "four fermions out of order": fc.Cdag[3]*fc.Cdag[1]*fc.C[2]*fc.C[0],
        "density times hopping": fc.N[2]*fc.Cdag[0]*fc.C[1],
        }
print()
for (label,op) in probes.items():
    d = np.max(np.abs(np.array(MO2matrix(op,obj)-MO2matrix(op.simplify(),obj))))
    print("%-28s terms %d -> %d, |matrix difference| = %.2e"%(label,
        len(op.op),len(op.simplify().op),d))
    assert d<1e-10
# and the anticommutator of two annihilation operators on different
# sites, which is zero only if the reordering sign is right
print("C[0]C[2] + C[2]C[0] proven zero:",
        (fc.C[0]*fc.C[2]+fc.C[2]*fc.C[0]).is_zero())

fig = plt.figure(figsize=(6.,4.2))
plt.plot(ns,np.array(tnum)*1e3,marker="o",label="numerical witness")
plt.plot(ns,np.array(tsym)*1e3,marker="s",label="canonical-form proof")
plt.yscale("log")
plt.xlabel("number of sites")
plt.ylabel("time to settle Hermiticity (ms)")
plt.title("Heisenberg chain, H = sum_i S_i . S_{i+1}")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("operator_canonical_form.png",dpi=150)
plt.show()
