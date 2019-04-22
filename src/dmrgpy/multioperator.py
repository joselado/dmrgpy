
import collections

# this class allows to define operators of the form
# A_0@A_1@....
# for mpscpp.x


class MultiOperator():
    """
    Object to deal with multioperators in mpscpp.x
    """
    def __init__(self,name="multioperator",c=1.0,term=True): # do nothing
        self.op = [] # empty list of sums of products
        self.name = name
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
        if type(a)!=MultiOperator: raise # not implemented
        out = self.copy() # create a copy
        out.op = self.op + a.op # sum the two operators
        out.i = self.i + a.i + 1 # increase the index
        out.clean()
        return out # return the sum
    def __rmul__(self,a):
        """Multiply by a number"""
        out = self.copy()
        for i in range(len(out.op)):
            out.op[i][0] = out.op[i][0]*a # multiply
        return out
    def __mul__(self,a):
        """Compute the product between two multioperators"""
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




