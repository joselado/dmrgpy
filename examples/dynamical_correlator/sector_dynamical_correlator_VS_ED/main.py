# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Regression test + visual check for submode="SECTOR", the sector-resolved
# Lehmann sum (src/dmrgpy/sectordc.py).
#
# The method sums <0|A|n><n|B|0> over the lowest eigenstates of the ONE
# quantum-number sector B|gs> lands in -- here B=S+, so the intermediate
# states are the Sz=+2 (2*Sz units, i.e. Sz=+1 physically) eigenstates and
# nothing else. Two things this example is meant to show:
#
#   1. Where the truncated Lehmann sum is complete, it is not an
#      approximation at all: it reproduces the exact ED correlator to
#      machine-ish precision, unlike KPM, which reconstructs the same
#      spectrum from a finite Chebyshev expansion and gets the peak
#      positions right but not the lineshape. That is the assert below.
#   2. What "complete" means is measurable rather than hoped for. The
#      right panel sweeps nex and plots the captured spectral weight,
#      sum_n |<n|B|0>|^2 / <0|B^dag B|0>, which is the method's own
#      statement of how much of the response it is describing.
#
# Both backends with sector support are checked (itensor_version=3 and
# "python"); no other backend has quantum numbers at all.
import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain
from dmrgpy.operatornames import name2MO

n = 6
es = np.linspace(-0.2,4.0,400)
delta = 0.1

def make_chain(backend):
    sc = spinchain.Spin_Chain([2 for i in range(n)])
    if backend=="python": sc.setup_python()
    else: sc.setup_cpp(backend)
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    sc.maxm = 40 ; sc.nsweeps = 8
    return sc,h

# --- exact reference: the full Lehmann sum over the whole spectrum ---
ref,h = make_chain("python")
Sp,Sm = name2MO("Sp",ref),name2MO("Sm",ref)
name = (Sm[0],Sp[0])   # <S-_0 ... S+_0>, the transverse channel
x_ed,y_ed = ref.get_dynamical_correlator(mode="ED",submode="ED",name=name,
                                          es=es,delta=delta)
y_ed = np.array(y_ed).real

spectra = {}
for backend in [3,"python"]:
    sc,_ = make_chain(backend)
    x,y,info = sc.get_dynamical_correlator(submode="SECTOR",name=name,
            nex=100,es=es,delta=delta,return_poles=True)
    spectra[backend] = (x,np.array(y).real,info)
    print("backend %-6s target sector %s, %d states, captured %.6f"
          %(str(backend),info["target_sector"],info["nex"],info["captured"]))
    diff = np.max(np.abs(np.array(y).real-y_ed))
    print("   max|SECTOR-ED| = %.2e (spectrum scale %.3f)"%(diff,np.max(np.abs(y_ed))))
    assert info["captured"]>1-1e-6, \
        "the sum should be complete here (nex covers the whole sector)"
    assert diff<1e-3*max(1.0,np.max(np.abs(y_ed))), \
        "backend %s disagrees with the exact ED correlator by %g"%(backend,diff)

# --- KPM on the same chain, for contrast: same physics, broadened ---
sc_kpm,_ = make_chain("python")
x_kpm,y_kpm = sc_kpm.get_dynamical_correlator(submode="KPM",name=name,
                                               es=es,delta=delta)

# --- how the captured weight (and the spectrum) fills in with nex ---
sc_conv,_ = make_chain("python")
nexs = [1,2,3,5,8,15]
captured,curves = [],{}
for nex in nexs:
    x,y,info = sc_conv.get_dynamical_correlator(submode="SECTOR",name=name,
            nex=nex,es=es,delta=delta,return_poles=True,quiet=True)
    captured.append(info["captured"])
    curves[nex] = np.array(y).real
    print("nex=%3d -> captured %.6f"%(nex,info["captured"]))
assert captured[-1]>captured[0], "captured weight should grow with nex"
print("TEST PASSED")

# --- plots ---
fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,4))
fig.subplots_adjust(wspace=0.3,bottom=0.15)

ax1.plot(x_ed,y_ed,c="black",lw=3,alpha=0.4,label="ED (exact Lehmann)")
for backend,c in ((3,"blue"),("python","green")):
    x,y,_ = spectra[backend]
    ax1.plot(x,y,c=c,ls="--",label="SECTOR, itensor_version=%s"%str(backend))
ax1.plot(x_kpm,np.array(y_kpm).real,c="red",ls=":",label="KPM")
ax1.set_xlabel("frequency [J]") ; ax1.set_ylabel(r"$S^{+-}(\omega)$")
ax1.set_title("Sz=+1 channel, n=%d Heisenberg"%n) ; ax1.legend(fontsize=8)

for nex in nexs:
    ax2.plot(es,curves[nex],label="nex=%d"%nex)
ax2.plot(x_ed,y_ed,c="black",lw=2,alpha=0.4,label="ED")
ax2.set_xlabel("frequency [J]") ; ax2.set_ylabel(r"$S^{+-}(\omega)$")
ax2.set_title("filling in with more states") ; ax2.legend(fontsize=7)

ax3.plot(nexs,captured,marker="o",c="purple")
ax3.axhline(1.0,c="black",ls="--",lw=1)
ax3.set_xlabel("nex (states kept in the sector)")
ax3.set_ylabel(r"captured weight $\sum_n |\langle n|B|0\rangle|^2/\langle B^\dagger B\rangle$")
ax3.set_title("the method's own error estimate") ; ax3.set_ylim(0,1.05)

plt.savefig("sector_dynamical_correlator_vs_ed.png",dpi=150)
print("Plot saved to sector_dynamical_correlator_vs_ed.png")
plt.show()
