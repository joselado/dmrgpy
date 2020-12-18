import numbers
import types
import collections
import numpy as np

# this class allows to define operators of the form
# A_0@A_1@....
# for mpscpp.x

def isnumber(x):
#    print(type(x),isinstance(x, numbers.Number))
    return isinstance(x, numbers.Number) or np.iscomplex(x)


ampo_counter = 0
use_jordan_wigner = True


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
    def simplify(self):
        from .multioperatortk import sympymultioperator
        return sympymultioperator.simplifyMO(self)
    def get_bandwidth(self,MBO):
        """Get the bandwidth"""
        return MBO.bandwidth(self)
    def is_hermitian(self):
        """Check if an operator is hermitian"""
        dh = self - self.get_dagger() ; dh = dh.simplify()
        return dh==0
    def is_antihermitian(self):
        """Check if an operator is antiHermitian"""
        dh = self + self.get_dagger() ; dh = dh.simplify()
        return dh==0
    def is_zero(self):
        return self.simplify()==0
    def copy(self):
        from copy import deepcopy
        return deepcopy(self) # return a copy
    def get_dagger(self):
        return get_dagger(self)
    def __neg__(self):
        return (-1)*self
    def __sub__(self,a):
        return self + (-1)*a
    def __rsub__(self,a):
        return -self + a
    def __add__(self,a):
        """Sum operation"""
        if a is None: return self.copy() # return the Hamiltonian
        elif type(a)==MultiOperator: # if it is a multioperator
          out = self.copy() # create a copy
          out.op = self.op + a.op # sum the two operators
          out.i = self.i + a.i + 1 # increase the index
          out.clean()
          return out # return the sum
        elif isnumber(a): # if it is a number
            return self+a*identity() # return identity
        else: 
            print(type(a),a)
            raise
    def __radd__(self,a): return self.__add__(a)
    def __rmul__(self,a):
        """Multiply by a number"""
        return self.__mul__(a)
    def __truediv__(self,a):
        if isnumber(a): # number
          return (1./a)*self
        else: raise
    def multiply_scalar(self,a):
        out = self.copy()
        for i in range(len(out.op)):
            out.op[i][0] = out.op[i][0]*a # multiply
        return out
    def __mul__(self,a):
        """Compute the product between two multioperators"""
        if type(a)==MultiOperator: return self.multiply_MO(a)
        elif type(a)==np.ndarray: raise  # prevent using rmul in array
        elif isnumber(a): return self.multiply_scalar(a)
        else: return NotImplemented
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
        if use_jordan_wigner: m = jordan_wigner(self)
        else: m = self
        write(m,name)
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


def identity():
    op = MultiOperator(term=True,c=1.0)
    op.add_operator("Id",1)
    return op

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
    write_ampo(MO2list(MO),name) # write in a file


def write_ampo(out,name):
    f = open(name,"w")
    f.write(str(len(out))+"\n") # number of lines
    for o in out:
      n = (len(o)-2)//2 # number of terms
      if n>=200: raise # C++ code needs to be recompiled
      f.write(str((len(o)-2)//2)+"\n") # number of terms
      for io in o:
          f.write(str(io)+"  ")
      f.write("\n")
    f.close()



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


def jordan_wigner(MO):
    """Use Jordan Wigner transformationin a multioperator"""
    from .multioperatortk import jordanwigner
    m = 0 # initialize output
    for ii in range(len(MO.op)): # loop over operators
        opi = MO.op[ii] # take this term
        n = len(opi)-1 # number of terms in the product
        c = opi[0] # take the complex value
        mi = list2MO(opi) # default output
        ls = [opi[kk+1][0] for kk in range(n)] # names of operators
        try:
          if (n==1): # one point operator
              i = opi[1][1]
              if "C" in ls or "Cdag" in ls:
                  mi = c*jordanwigner.one_fermion(ls[0],i)
          elif (n==2): # two point operators
              i,j = opi[1][1],opi[2][1]
              if "C" in ls or "Cdag" in ls:
                  mi = c*jordanwigner.two_fermions(ls[0],i,ls[1],j)
          elif (n==4): # four point operators
              i,j,k,l = opi[1][1],opi[2][1],opi[3][1],opi[4][1]
              if "C" in ls or "Cdag" in ls:
                  mi = c*jordanwigner.four_fermions(ls[0],i,ls[1],j,ls[2],
                          k,ls[3],l) 
          else: 
              if "C" in ls or "Cdag" in ls: raise # not implemented
        except: # brute force 
            inds = [opi[kk+1][1] for kk in range(n)]
            mi = c*jordanwigner.product2jw(ls,inds)
        m = m + mi
    return m



def get_dagger(self,conjugate=True):
    """Return the dagger of a multioperator"""
    m = 0 # initialize
    for opi in self.op: # loop over terms
        n = len(opi) # number of terms in the product
        c = opi[0] # coefficient
        mi = np.conjugate(c) # initialize
        for i in range(1,n): # loop over terms in the product
            name = opi[i][0] # name
            jj = opi[i][1] # index
            if name=="C": name2="Cdag"
            elif name=="Cdag": name2="C"
            elif name=="A": name2="Adag"
            elif name=="Adag": name2="A"
            elif name=="Sp": name2="Sm"
            elif name=="S+": name2="S-"
            elif name=="Sm": name2="Sp"
            elif name=="S-": name2="S+"
            elif name=="Sig": name2="SigDag"
            elif name=="Tau": name2="TauDag"
            elif name=="SigDag": name2="Sig"
            elif name=="TauDag": name2="Tau"
            else: name2 = name
            mi = obj2MO([[name2,jj]])*mi
        m = m + mi # add contribution
    return m # return MO








