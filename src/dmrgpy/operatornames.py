



def recognize(name):
  """
  Recognize which correlator you want to compute
  """
  # single operators, for expectation values
  if name=="X": return "Id","Sx"
  if name=="Y": return "Id","Sy"
  if name=="Z": return "Id","Sz"
  # next cases
  spins = [] # names for spin correlators
  for s1 in ["X","Y","Z","+","-"]:
    for s2 in ["X","Y","Z","+","-"]: spins += [s1+s2] # create all the possibilities
  # keep in mind that the KPM will use the dagger of the first operator!
  if name in spins: # spin correlator
    if name[0]=="X": namei="Sx"
    elif name[0]=="Y": namei="Sy"
    elif name[0]=="Z": namei="Sz"
    elif name[0]=="+": namei="Sp" # the other way around!
    elif name[0]=="-": namei="Sm" # the other way around!
    else: raise
    if name[1]=="X": namej="Sx"
    elif name[1]=="Y": namej="Sy"
    elif name[1]=="Z": namej="Sz"
    elif name[1]=="-": namej="Sm"
    elif name[1]=="+": namej="Sp"
    else: raise
#    if self.sites[i] !=1 or self.sites[j]!=1:
#        if name!="ZZ": raise  # fermions only accept ZZ
  else: # fermionic correlator
#    if self.sites[i] !=1 or self.sites[j]!=1: raise # only for fermions
    if name=="cdc": namei = "Cdag" ; namej = "C"
    elif name=="cdcup": namei = "Cdagup" ; namej = "Cup"
    elif name=="cdcdn": namei = "Cdagdn" ; namej = "Cdn"
    elif name=="ccd": namei = "C" ; namej = "Cdag"
    elif name=="cc": namei = "C" ; namej = "C"
    elif name=="deltadelta" or name=="delta":
        namei = "delta" ; namej = "delta"
    elif name=="deltadeltad":
        namei = "delta" ; namej = "deltad"
    elif name=="densitydensity":
        namei = "N" ; namej = "N"
    elif name=="density": # density density correlator
        namei = "N" ; namej = "N"
    else: raise
#  else:
#      print("Dynamical correlator not recognised")
#      raise
  return namei,namej


def hermitian(name):
    """Return the Hermitian name"""
    if name=="Sp": return "Sm"
    elif name=="Sm": return "Sp"
    elif name=="C": return "Cdag"
    elif name=="Cup": return "Cdagup"
    elif name=="Cdn": return "Cdagdn"
    elif name=="Cdagdn": return "Cdn"
    elif name=="Cdagup": return "Cup"
    elif name=="Cdag": return "C"
    return name


def name2MO(name,self):
    if name=="C": return self.C
    elif name=="Cdag": return self.Cdag
    elif name=="N": return self.N
    elif name=="Adag": return self.Adag
    elif name=="A": return self.A
    elif name=="Sx": return self.Sx
    elif name=="Sy": return self.Sy
    elif name=="Sz": return self.Sz
    # Sp/Sm have no native backend operator name (unlike Sx/Sy/Sz), so build
    # them from Sx/Sy on the fly: S+ = Sx + i*Sy, S- = Sx - i*Sy (matching
    # the Sp<->Sm dagger convention in multioperator.py/sympymultioperator.py)
    elif name=="Sp": return [self.Sx[i]+1j*self.Sy[i] for i in range(self.ns)]
    elif name=="Sm": return [self.Sx[i]-1j*self.Sy[i] for i in range(self.ns)]
    else:
        raise ValueError("Unrecognized operator name "+str(name))

def str2MO(self,name,i=0,j=0):
    """Normalize the `name=` argument of a correlator into a pair of
    MultiOperators: either a documented string ("ZZ", "cdc", ...) plus
    site indices i/j, or an explicit (A,B) MultiOperator pair.

    The string branch used to cover Sx/Sy/Sz/Sp/Sm only and end in a bare
    `raise`, so the documented fermionic names recognize() itself returns
    ("cdc", "ccd", "densitydensity", ...) died with the opaque
    "RuntimeError: No active exception to reraise" -- including from
    fermionchain.get_gr, which passes name="cdc" itself. name2MO covers
    the remaining families, so the string form now works wherever the
    operator exists on the chain, and says what is wrong when it does
    not."""
    from . import multioperator
    if type(name)==str:
        n1,n2 = recognize(name)
        def f(n,k):
            if n=="Id": return multioperator.identity()
            try: ops = name2MO(n,self)
            except (ValueError,AttributeError):
                raise ValueError(
                    "Operator "+repr(n)+" (from correlator name "+repr(name)
                    +") is not available on this chain; pass an explicit "
                    "(A,B) MultiOperator pair instead")
            return ops[k]
        return [f(n1,i),f(n2,j)]
    elif type(name[0])==multioperator.MultiOperator and type(name[1])==multioperator.MultiOperator:
        return [name[0],name[1]]
    else:
        raise ValueError(
            "name must be a documented correlator string (with i=/j=) or a "
            "pair of MultiOperators, got "+repr(name))

