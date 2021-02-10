from sympy import *
import numbers
from .. import multioperator as classicMO

class MultiOperator(Symbol):
    """Multioperator class"""
    def __init__(self,name=None,c=1.0,O="Id_1"):
        super().__init__(O,commutative=False) # initialize
        if name is None:
            self.name = "ampo_operator_"+str(ampo_counter)
            ampo_counter += 1
        else: self.name = name
        return self*c
    def write(self,*args):
        """Write in file"""
        return classicMO.write_ampo(symbol2MO(self,*args))


def simplifyMO(self):
    a = MO2symbol(self)
    a = simplify(a)
    return symbol2MO(a)



def MO2symbol(self):
    """Transform the custom MO to a sympy symbol"""
    out = 0
    for o in self.op:
        tmp = o[0] # numerical value
        for i in range(1,len(o)):
            name = str(o[i][0])+"_"+str(o[i][1])
            tmp = tmp*Symbol(name,commutative=False)
        out = out + tmp
    return out

def symbol2MO(t):
    """Transform a sympy symbol to a MO"""
    from ..multioperator import list2MO
    if type(t)==Add: # right type
        out = 0
        for o in t.args:
            out = out + symbol2MO(o)
        return out
    elif type(t)==Symbol: 
        name = t.name.split("_")
        return list2MO([1.0,[name[0],int(name[1])]]) # return the MO
    elif type(t)==Mul: # multiplication
        out = 1
        for ti in t.args:
            out = out*symbol2MO(ti) # multiply
        return out
    elif type(t)==Pow: # power
        out = 1
        if len(t.args)>2: raise # unrecognized
        for i in range(t.args[1]): # loop
            out = out*symbol2MO(t.args[0]) # multiply
        return out
    elif type(t)==Integer: return complex(t)
#    elif type(t)==numbers.ImaginaryUnit: return 1j
    elif type(t)==Float: return complex(t)
    elif hasattr(t, 'is_number'): # if it is a sumpy number
        if t.is_number: return complex(t)
    else: raise




dagdict = dict()
dagdict["C"] = "Cdag"
dagdict["Cdag"] = "C"
dagdict["A"] = "Adag"
dagdict["Adag"] ="A"
dagdict["Sp"] = "Sm"
dagdict["S+"] = "S-"
dagdict["Sm"] = "Sp"
dagdict["S-"] = "S+"
dagdict["Sig"] = "SigDag"
dagdict["Tau"] = "TauDag"
dagdict["SigDag"] = "Sig"
dagdict["TauDag"] = "Tau"

def get_dagger(self):
    return self.copy()
#    for d in dagdict

