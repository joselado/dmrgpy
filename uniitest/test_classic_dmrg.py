import unittest
import numpy as np
from dmrgpy import simplechains

error = 1e-2 # acceptable accuracy

class TestDMRG(unittest.TestCase):

    def test_gs_dmrg(self):
        # compute energy with MPS and classic DMRG
        b = np.random.random(3) 
        ssc = simplechains.SSC(s=.5,n=10,b=b)
        de = ssc.gs_energy(mode="MPS") - ssc.gs_energy(mode="classicDMRG")
        self.assertTrue(abs(de)<error)
    def test_excited_dmrg(self):
        # compute excitation energies
        b = np.random.random(3) ; ne = 4
        ssc = simplechains.SSC(s=.5,n=10,b=b)
        e0 = ssc.get_excited(mode="MPS",n=ne) 
        e1 = ssc.get_excited(mode="classicDMRG",n=ne)
        de = e0 - e1
        de = np.sum(np.abs(de))
        passed = de < error
        if not passed:
            print("ERROR in energies is",de)
            print("Energies with MPS",e0)
            print("Energies with classic DMRG",e1)
        self.assertTrue(passed)

if __name__ == '__main__':
    unittest.main()
