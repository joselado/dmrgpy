



def recognize(self,name):
  fcorr = ["cdc","cdcup","ccd","cdcdn","deltadelta","deltadeltad"] # correlators for fermions
  spins = [] # names for spin correlators
  for s1 in ["X","Y","Z","+","-"]:
    for s2 in ["X","Y","Z","+","-"]: spins += [s1+s2] # create all the possibilities
  # keep in mind that the KPM will use the dagger of the first operator!
  if name in spins: # spin correlator
    if name[0]=="X": namei="Sx"
    elif name[0]=="Y": namei="Sy"
    elif name[0]=="Z": namei="Sz"
    elif name[0]=="+": namei="Sm" # the other way around!
    elif name[0]=="-": namei="Sp" # the other way around!
    else: raise
    if name[1]=="X": namej="Sx"
    elif name[1]=="Y": namej="Sy"
    elif name[1]=="Z": namej="Sz"
    elif name[1]=="-": namej="Sm"
    elif name[1]=="+": namej="Sp"
    else: raise
#    if self.sites[i] !=1 or self.sites[j]!=1:
#        if name!="ZZ": raise  # fermions only accept ZZ
  elif name in fcorr: # fermionic correlator
    if self.sites[i] !=1 or self.sites[j]!=1: raise # only for fermions
    if name=="cdc": namei = "Cdag" ; namej = "Cdag"
    elif name=="cdcup": namei = "Cdagup" ; namej = "Cdagup"
    elif name=="cdcdn": namei = "Cdagdn" ; namej = "Cdagdn"
    elif name=="ccd": namei = "C" ; namej = "C"
    elif name=="deltadelta" or name=="delta":
        namei = "delta" ; namej = "delta"
    elif name=="deltadeltad":
        namei = "delta" ; namej = "deltad"
    else: raise
  else:
      print("Dynamical correlator not recognised")
      raise
  return namei,namej

