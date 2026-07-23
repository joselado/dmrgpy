import numpy as np
from .. import multioperator
from ..algebra import algebra

# Full-spectrum ED building blocks for the third-order Kondo/STM
# perturbation theory (Ternes, arXiv:1505.04430). Unlike the ground-state
# ED already used elsewhere in dmrgpy (EDchain.get_gs_array()), the third
# order Kondo term sums over *every* eigenstate as a virtual intermediate
# state, so this always needs the full spectrum, not just the low-energy
# states DMRG or a Lanczos-type ED would normally target.


class KondoSpectrum():
    """Holds everything derived once from a chain's full ED spectrum that
    the second- and third-order conductance terms need: eigenenergies
    (shifted so the ground state is at 0), Boltzmann occupations at
    temperature T, and the impurity spin operators Sx/Sy/Sz transformed
    into the Hamiltonian eigenbasis."""
    def __init__(self, chain, site, T, kB=8.617333262e-5):
        # kB in eV/K by default (so T is given in Kelvin, energies in eV,
        # matching how the rest of dmrgpy's examples quote energies); pass
        # kB=1 if you'd rather work in units where T is already an energy
        if T <= 0.: raise ValueError("KondoSpectrum requires T>0")
        self.chain = chain
        self.site = site
        self.T = T
        self.kB = kB
        edobj = chain.get_ED_obj() # full-spectrum ED backend (pychain)
        emu, vs = edobj.get_diagonalized_hamiltonian() # all eigenpairs
        emu = np.array(emu, dtype=float)
        order = np.argsort(emu) # eigh output should already be sorted
        emu = emu[order]
        vs = np.array(vs)[:, order]
        self.e = emu - emu[0] # ground state energy set to 0
        self.vs = vs
        self.dim = len(self.e)
        kT = kB*T
        p = np.exp(-self.e/kT)
        self.p = p/np.sum(p) # Boltzmann occupations, eq. "boltzmann"
        self.Sx = self._eigenbasis_spin_operator("Sx")
        self.Sy = self._eigenbasis_spin_operator("Sy")
        self.Sz = self._eigenbasis_spin_operator("Sz")
    def _eigenbasis_spin_operator(self, name):
        edobj = self.chain.get_ED_obj()
        op = getattr(self.chain, name)[self.site] # MultiOperator, this site
        mat = algebra.todense(multioperator.MO2matrix(op, edobj))
        mat = np.array(mat, dtype=complex)
        return np.conjugate(self.vs.T)@mat@self.vs # Xi^alpha_{ab} = <a|S|b>
