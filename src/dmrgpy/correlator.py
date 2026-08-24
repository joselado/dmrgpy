from . import operatornames
import numpy as np


def get_correlator(self,pairs=[[]],name="SS",**kwargs):
    """Compute a certain static correlator in the spin chain"""
    if name=="SS":
        def getop(i,j): 
            m = self.Sx[i]*self.Sx[j]
            m = m + self.Sy[i]*self.Sy[j]
            m = m + self.Sz[i]*self.Sz[j]
            return m
    else:
      namei,namej = operatornames.recognize(name) # return that one
      opi = operatornames.name2MO(namei,self)
      opj = operatornames.name2MO(namej,self)
      def getop(i,j): 
          return opi[i]*opj[j]
    return np.array([self.vev(getop(i,j),**kwargs) for (i,j) in pairs])
