
import collections
import numpy as np

# this class allows to define operators of the form
# A_0@A_1@....
# for mpscpp.x


ampo_counter = 0



class MultiOperator():
    """
    Object to deal with multioperators in mpscpp.x
    """
    def __init__(self,name=None,c=1.0,term=True): # do nothing
        global ampo_counter
        self.op = [] # empty list of sums of products
        if name is None:
            self.name = "ampo_operator_"+str(ampo_counter)
            ampo_counter += 1
        else: self.name = name
        self.i = -1 # initialize
        if term: self.new_term(c=c) # generate the first term
    def add_operator(self,name,i):
        """Store operator"""
        self.op[self.i].append([name,i]) # append that name
    def new_term(self,c=1.0):
        """Add a new term"""
        self.i += 1 # increase the counter
        self.op.append([c]) # initialize
    def copy(self):
        from copy import deepcopy
        return deepcopy(self) # return a copy
    def __add__(self,a):
        """Sum operation"""
        if a is None: return self.copy() # return the Hamiltonian
        if not type(a)==MultiOperator: # assume that it is a number
            if np.abs(a)==0.0: return self.copy() # return object
            else: raise
        else: # multioperator
          out = self.copy() # create a copy
          out.op = self.op + a.op # sum the two operators
          out.i = self.i + a.i + 1 # increase the index
          out.clean()
          return out # return the sum
    def __radd__(self,a): return self.__add__(a)
    def __rmul__(self,a):
        """Multiply by a number"""
        return self.__mul__(a)
    def multiply_scalar(self,a):
        out = self.copy()
        for i in range(len(out.op)):
            out.op[i][0] = out.op[i][0]*a # multiply
        return out
    def __mul__(self,a):
        """Compute the product between two multioperators"""
        if type(a)==MultiOperator: return self.multiply_MO(a)
        else: return self.multiply_scalar(a)
    def multiply_MO(self,a):
        out = self.copy() # copy operator
        out.i = (self.i+1)*(a.i+1) # total number of terms
        out.op = [] # empty list
        for io in self.op: # loop over first operator
          for jo in a.op: # loop over second operator
              o = [io[0]*jo[0]] # compute coefficient
              o = o + [io[i] for i in range(1,len(io))] # add
              o = o + [jo[i] for i in range(1,len(jo))] # add
              out.op.append(o) # store contribution
        out.clean()
        return out # return operator
    def clean(self):
        """Remove terms with zero weight"""
        op = []
        for o in self.op:
            if abs(o[0])>1e-8: op.append(o) # store
        self.i = self.i - (len(self.op)-len(op)) # redefine
        self.op = op # redefine
    def write(self,name=None):
        """Write in a file"""
        if name is None: name = self.name+".in"
        write(self,name)
    def get_dict(self):
        """Return the dictionary to be used in tasks.in"""
        d = dict()
        ii = 0
        d[self.name+"_n"] = str(len(self.op)) # write the coefficient
        for iop in self.op: # loop over operators
          name = self.name+"_operator_"+str(ii)
          d[name+"_coefficient_real"] = str(iop[0].real) # write the coefficient
          d[name+"_coefficient_imag"] = str(iop[0].imag) # write the coefficient
          d[name+"_nterms"] = len(iop)-1 # write the coefficient
          ii += 1 # increase counter
          for i in range(len(iop)-1): # loop over terms of the product
              o = iop[i+1] # get the term
              name0 = name+"_term_"+str(i)+"_name"
              name1 = name+"_term_"+str(i)+"_site"
              d[name0] = o[0]
              d[name1] = o[1]
        return d


def zero():
    return MultiOperator(term=True,c=0.0)


def MO2list(self):
    """Conver a multioperator into a list"""
    out = []
    for iop in self.op: # loop over operators
        o = []
        o.append(iop[0].real) # real part
        o.append(iop[0].imag) # imaginary part
        for i in range(len(iop)-1): # loop over terms
          otmp = iop[i+1]
          o.append(otmp[0])
          o.append(otmp[1]+1)
        out.append(o)
    return out


def write(MO,name):
    """Write a multioperator in a file"""
    from .ampotk import write_ampo
    write_ampo(MO2list(MO),name) # write in a file



def obj2MO(a,name="multioperator"):
    """
    Convert an input in a multioperator
    """
    if isinstance(a, collections.Iterable): # if it is a tuple
        mo = MultiOperator(name=name) # create object
        for ia in a:
            mo.add_operator(ia[0],ia[1])
        return mo
    elif type(a)==MultiOperator:
        a.name = name
        return a
    else: raise # unidentified input



def list2MO(l):
    """Convert a list into a multioperator"""
    mo = MultiOperator(c=l[0]) # create output
    for i in range(1,len(l)):
        mo.add_operator(l[i][0],l[i][1])
    return mo




def MO2vijkl(mo):
    """Extract from a multioperator the part corresponding to
    Vijkl interaction"""
    mout = MultiOperator(term=False) # create output
    d = dict() # dictionary
    for o in mo.op: # loop over terms in the sum
        if len(o)==5: # four products plus scalar
            if o[1][0]=="Cdag" and o[2][0]=="C" and o[3][0]=="Cdag" and o[4][0]=="C":
                ind = (o[1][1],o[2][1],o[3][1],o[4][1])
                d[ind] = o[0] # set the value of the coupling
        else: mout = mout + list2MO(o) # convert into a multioperator
    def f(i,j,k,l):
        try: return d[(i,j,k,l)]
        except: return 0.0
    if len(mout.op)==0: mout = None
    return f,mout



def MO2matrix(MO,obj):
    """Given a certain object containing the method "get_operator",
    return a matrix"""
    out = 0.0
    for iop in MO.op: # loop over components
        otmp = iop[0]*obj.get_identity() # factor
        for i in range(len(iop)-1): # loop over terms in the product
            term = iop[i+1] # get this term
            otmp = otmp@obj.get_operator(term[0],term[1]) # multiply
        out = out + otmp
    return out # return matrix






