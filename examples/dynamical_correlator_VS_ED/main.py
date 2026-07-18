# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

# Regression test: dynamical correlator <Sz_0(t) Sz_0(0)> submodes must
# agree with the exact ED reference (mode="ED", submode="ED" -- a plain
# Lehmann-representation sum over the exact spectrum, independently
# cross-checked here against submode="INV", the resolvent/matrix-
# inversion route: the two agree to ~1e-15, confirming both are valid
# ground truth). No existing dynamical_correlator example checks against
# ED with an assertion --
# examples/dynamical_correlator/dynamical_correlator_fermionic_ED_VS_DMRG
# and .../dynamical_correlator_fermionic_ED_VS_KPM_VS_CVM only plot (and
# the latter has a bare "raise # this should be double checked" left in
# it, i.e. it doesn't even run).
#
# Two submodes are checked very differently, because they are different
# kinds of approximation:
#
#   - CVM (correction-vector method, cvm.py -- rewritten as a pure-Python
#     CG solve after a prior version had real blow-ups, see project
#     memory) solves a resolvent linear system at each frequency, which
#     for an exact ground state is numerically exact (like ED's own
#     "INV" submode) rather than a controlled approximation -- confirmed
#     directly, CVM matches ED pointwise to ~1e-13 here. So CVM gets a
#     tight, pointwise tolerance.
#
#   - KPM (Chebyshev moment expansion, kpmdmrg.py) is a genuine
#     controlled approximation: it reconstructs the spectral function
#     from a finite number of moments convolved with a broadening
#     kernel, which does NOT reproduce the same lineshape as ED's fixed-
#     delta Lorentzian broadening point-by-point -- confirmed directly,
#     the pointwise difference stayed large even as the moment count
#     (kpm_n_scale) was raised well past its default. What DOES match
#     closely is where the physical excitation energies show up: the
#     peak positions of the reconstructed spectrum. So KPM is checked by
#     locating peaks in both spectra (scipy.signal.find_peaks) and
#     matching each ED peak to the nearest KPM peak within a tolerance,
#     rather than comparing values point-by-point.
import numpy as np
from scipy.signal import find_peaks
from dmrgpy import spinchain

n = 4
es = np.linspace(-0.5,4.0,200) # fine enough to resolve peak positions
delta = 0.15

def make_chain(itensor_version=None):
    sc = spinchain.Spin_Chain([2 for i in range(n)])
    if itensor_version is not None:
        if itensor_version!="python": sc.setup_cpp(itensor_version)
        else: sc.setup_python()
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h)
    return sc

def peak_energies(x,y):
    """Energies of the resolvable peaks (>=5% of the tallest one)"""
    y = np.array(y).real
    x = np.array(x)
    idx,_ = find_peaks(y,height=0.05*np.max(y))
    return x[idx]

# --- exact ED ground truth, cross-checked against itself ---
sc_ed = make_chain()
name = (sc_ed.Sz[0],sc_ed.Sz[0]) # <Sz_0(t) Sz_0(0)>
x_ed,y_ed = sc_ed.get_dynamical_correlator(mode="ED",submode="ED",name=name,es=es,delta=delta)
x_inv,y_inv = sc_ed.get_dynamical_correlator(mode="ED",submode="INV",name=name,es=es,delta=delta)
y_ed,y_inv = np.array(y_ed),np.array(y_inv)
diff_ed_inv = np.max(np.abs(y_ed-y_inv))
print("ED (Lehmann) vs ED (INV) max|diff| =",diff_ed_inv)
assert diff_ed_inv<1e-6, "the two exact ED submodes disagree by %g -- one of them is broken"%diff_ed_inv

peaks_ed = peak_energies(x_ed,y_ed)
print("ED peak energies =",peaks_ed)
assert len(peaks_ed)>=2, "expected at least 2 resolvable peaks in this window, found %d -- adjust es/delta"%len(peaks_ed)

# --- CVM: resolvent method, should be pointwise-exact ---
# tol chosen to actually discriminate: confirmed directly that the prior
# CVM implementation (hand-rolled BiCGSTAB in chain_session.h, replaced
# for real blow-ups on larger systems -- see project memory) still
# disagrees with ED by ~1e-6..1.6e-6 even on this small, well-conditioned
# system, while the current CG-based solver matches to ~1e-13..1e-14; a
# looser tolerance like 1e-3 would silently miss that regression.
tol_cvm = 1e-6
for v in [2,3,"python"]:
    sc = make_chain(v)
    name = (sc.Sz[0],sc.Sz[0])
    x_cvm,y_cvm = sc.get_dynamical_correlator(mode="DMRG",submode="CVM",name=name,es=es,delta=delta)
    y_cvm = np.array(y_cvm)
    diff = np.max(np.abs(y_cvm-y_ed))
    print("CVM(v%s) vs ED max|diff| = %.2e"%(v,diff))
    assert diff<tol_cvm, "v%s CVM disagrees with ED by %g (tol=%g)"%(v,diff,tol_cvm)

# --- KPM: Chebyshev moment expansion, only peak positions (excitation
# energies) are expected to match, not the pointwise lineshape ---
tol_peak = 0.1 # generous relative to the frequency spacing (es step
                # ~0.023) and delta=0.15 broadening used above
for v in [2,3,"python"]:
    sc = make_chain(v)
    name = (sc.Sz[0],sc.Sz[0])
    x_kpm,y_kpm = sc.get_dynamical_correlator(mode="DMRG",submode="KPM",name=name,es=es,delta=delta)
    peaks_kpm = peak_energies(x_kpm,y_kpm)
    print("KPM(v%s) peak energies ="%v,peaks_kpm)
    assert len(peaks_kpm)>=len(peaks_ed), \
        "v%s KPM found %d peaks, expected at least the %d ED found"%(v,len(peaks_kpm),len(peaks_ed))
    # match each ED peak to its nearest KPM peak
    for e_ed in peaks_ed:
        closest = peaks_kpm[np.argmin(np.abs(peaks_kpm-e_ed))]
        diff = abs(closest-e_ed)
        print("  ED peak at %.4f -> nearest KPM peak at %.4f (diff=%.3f)"%(e_ed,closest,diff))
        assert diff<tol_peak, "v%s KPM has no peak within %g of the ED peak at %.4f (closest: %.4f)"%(v,tol_peak,e_ed,closest)

print("TEST PASSED")
