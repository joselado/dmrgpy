import numbers
import types
import collections
import numpy as np

# this class allows to define operators of the form
# A_0@A_1@....

# Concrete scalar types checked first (fast path): isinstance() against the
# numbers.Number ABC has to walk its __instancecheck__ machinery, which is
# measurably slower than a plain isinstance() against a tuple of concrete
# types, and this is called on every MultiOperator +/-/*. Uncommon numeric
# types (e.g. exotic numpy/Python Number subclasses) still fall back to the
# ABC check below.
_fast_number_types = (int, float, complex, np.integer, np.floating,
                       np.complexfloating)

def isnumber(x):
    return isinstance(x, _fast_number_types) or isinstance(x, numbers.Number)


ampo_counter = 0
use_jordan_wigner = True
n_mpo_max = 100 # maximum number of MPO products
clean_threshold = 1e-8 # coefficient magnitude below which a term is dropped


def _filter_small(op):
    """Return a new term list with near-zero-coefficient terms dropped,
    without mutating the input (terms may be shared between several
    MultiOperator objects, see MultiOperator.copy)."""
    return [o for o in op if abs(o[0])>clean_threshold]


class MultiOperator():
    """
    Object to deal with multioperators
    """
    __slots__ = ("op","name","i")
    def __init__(self,name=None,c=1.0,term=True): # do nothing
        global ampo_counter
        self.op = [] # empty list of sums of products
        if name is None:
            self.name = "ampo_operator_"+str(ampo_counter)
            ampo_counter += 1
        else: self.name = name
        # self.i could be eventually removed, as it does nothing
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
    def max_op_prod(self):
        """Maximum number of oprators in a product"""
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
        """Return a cheap copy of this operator.

        Terms are never mutated in place once they are part of self.op
        (multiply_scalar, multiply_MO, etc. always build fresh term
        lists instead), so sharing the term objects between the
        original and the copy and only duplicating the outer list is
        safe, and much cheaper than copy.deepcopy.
        """
        out = MultiOperator.__new__(MultiOperator)
        out.op = list(self.op)
        out.name = self.name
        out.i = self.i
        return out
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
          out.i = len(out.op)-1 # increase the index
          # No clean() here: concatenating two already-clean term lists
          # cannot introduce new near-zero terms, and calling clean()
          # (an O(len(op)) rescan) on every "+" is what made building a
          # Hamiltonian one term at a time (H = H + term, used pervasively
          # across this codebase) quadratic. Filtering happens lazily,
          # once, right before terms are consumed (write()/to_terms()).
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
        if not isnumber(a): raise # number
        out = self.copy()
        # Build fresh term lists rather than mutating out.op[i][0] in
        # place: out.copy() shares the term objects with self (see
        # copy()), so mutating them here would silently corrupt self too.
        out.op = [[o[0]*a]+o[1:] for o in self.op]
        return out
    def __mul__(self,a):
        """Compute the product between two multioperators"""
        if type(a)==MultiOperator: return self.multiply_MO(a)
        elif type(a)==np.ndarray: raise  # prevent using rmul in array
        elif isnumber(a): return self.multiply_scalar(a)
        else: return NotImplemented
    def multiply_MO(self,a):
        # self.op gets fully replaced below, so there is no point copying
        # it first (the old code's self.copy() here just deepcopied self.op
        # only to immediately discard it) - build a fresh, empty object.
        out = MultiOperator.__new__(MultiOperator)
        out.name = self.name
        sop,aop = self.op,a.op
        out.op = [[io[0]*jo[0]]+io[1:]+jo[1:] for io in sop for jo in aop]
        out.i = len(out.op)-1
        return out # return operator
    def clean(self):
        """Remove terms with zero weight"""
        self.op = _filter_small(self.op)
        self.i = len(self.op)-1
    def write(self,name=None):
        """Write in a file"""
        if name is None: name = self.name+".in"
        if use_jordan_wigner: m = jordan_wigner(self)
        else: m = self
        write(m,name)
    def to_terms(self):
        """
        Return this operator as [(coeff, [(name,1-based-site),...]), ...],
        for the in-process pybind11 extension (mpscpp2/bindings.cc). This is
        exactly what write()/get_dict() already compute (Jordan-Wigner
        transform, 1-based site indices), just returned as Python objects
        instead of being written to a file/tasks.in dict.
        """
        if use_jordan_wigner: m = jordan_wigner(self)
        else: m = self
        return [(complex(term[0]),[(name,i+1) for (name,i) in term[1:]])
                for term in _filter_small(m.op)]
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


def msum(ops):
    """Sum an iterable of MultiOperator objects (bare 0/None entries are
    treated as no-ops) in a single linear pass over all their terms.

    This is the bulk equivalent of writing "h = h + term" in a loop: doing
    that N times is O(N^2) (each "+" copies/concatenates the whole term
    list built so far), while collecting the terms and summing once here
    is O(total number of terms). Prefer this in any loop that assembles
    an operator out of many pieces (long-range couplings, sums over
    site/orbital indices, ...).
    """
    terms = []
    name = None
    for o in ops:
        if o is None: continue
        if isnumber(o):
            if o==0: continue
            o = o*identity()
        if name is None: name = o.name
        terms.extend(o.op)
    out = MultiOperator(term=False)
    if name is not None: out.name = name
    out.op = terms
    out.i = len(terms)-1
    return out


def identity():
    op = MultiOperator(term=True,c=1.0)
    op.add_operator("Id",1)
    return op

def MO2list(self):
    """Convert a multioperator into a list"""
    out = []
    for iop in _filter_small(self.op): # loop over operators, dropping
        o = []                         # any near-zero-weight ones
        o.append(iop[0].real) # real part
        o.append(iop[0].imag) # imaginary part
        for otmp in iop[1:]: # loop over terms
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
      if n>=n_mpo_max: 
          print("Too long MPO")
          raise # C++ code needs to be recompiled
      f.write(str((len(o)-2)//2)+"\n") # number of terms
      for io in o:
          f.write(str(io)+"  ")
      f.write("\n")
    f.close()



def obj2MO(a,name="multioperator"):
    """
    Convert an input in a multioperator
    """
    if isinstance(a, (list,tuple)): # if it is a list/tuple of [name,site]
        mo = MultiOperator(name=name) # create object
        for ia in a:
            mo.add_operator(ia[0],ia[1])
        return mo
    elif type(a)==MultiOperator:
        a.name = name
        return a
    elif isnumber(a):
        out = a*identity()
        out.name = name # a bare number has no name of its own to keep
        return out
    else: raise # unidentified input



def list2MO(l):
    """Convert a list into a multioperator"""
    mo = MultiOperator(c=l[0]) # create output
    for i in range(1,len(l)):
        mo.add_operator(l[i][0],l[i][1])
    return mo




def MO2matrix(MO,obj):
    """Given a certain object containing the method "get_operator",
    return a matrix"""
    out = 0.0*obj.get_identity() # initialize
    # This is a consumption point (like write()/to_terms()), so drop
    # near-zero-coefficient terms here too -- e.g. the "0*identity()"
    # placeholder that "0 + MultiOperator" (the "h = 0; h = h + term"
    # idiom used pervasively across this codebase) creates via __radd__.
    # Without this, a term naming an operator/site combination that was
    # never registered with obj (e.g. Parafermionic_Chain's ED backend
    # has no ("Id",1) entry) crashes get_operator() below even though
    # its contribution is exactly zero.
    for iop in _filter_small(MO.op): # loop over components
        otmp = iop[0]*obj.get_identity() # factor
        for i in range(len(iop)-1): # loop over terms in the product
            term = iop[i+1] # get this term
            otmp = otmp@obj.get_operator(term[0],term[1]) # multiply
        out = out + otmp
    return out # return matrix


_spinful_fermion_names = {"Cup","Cdagup","Cdn","Cdagdn"}


def jordan_wigner(MO):
    """Use Jordan Wigner transformationin a multioperator"""
    from .multioperatortk import jordanwigner
    from .multioperatortk import jordanwigner_spinful
    # Accumulate into a flat Python list and wrap once at the end, instead
    # of "m = m + mi" per term: the latter goes through the public "+"
    # (copy + list concat) once per term of MO.op, which turns this loop
    # into an accidental O(n^2) even though the work is inherently O(n).
    outterms = []
    for opi in MO.op: # loop over operators
        n = len(opi)-1 # number of terms in the product
        c = opi[0] # take the complex value
        mi = None # None means "term unchanged", set below if JW applies
        ls = [opi[kk+1][0] for kk in range(n)] # names of operators
        if any(l in _spinful_fermion_names for l in ls):
            # Native spinful-fermion sites (fermionchain.py's
            # Spinful_Fermionic_Chain_Native): always use the generic
            # per-factor dressing, there is no separate optimized 2-/
            # 4-point path for this operator family (see
            # jordanwigner_spinful.py's docstring for why the generic
            # recipe is already exact, not an approximation).
            inds = [opi[kk+1][1] for kk in range(n)]
            mi = c*jordanwigner_spinful.product2jw(ls,inds)
        else:
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
        if mi is None: outterms.append(list(opi)) # term unchanged
        else: outterms.extend(mi.op) # flatten in the transformed terms
    out = MultiOperator(term=False)
    out.op = outterms
    out.i = len(outterms)-1
    return out



_dagger_name = {"C":"Cdag","Cdag":"C","A":"Adag","Adag":"A",
                "Cup":"Cdagup","Cdagup":"Cup","Cdn":"Cdagdn","Cdagdn":"Cdn",
                "Aup":"Adagup","Adagup":"Aup","Adn":"Adagdn","Adagdn":"Adn",
                "Sp":"Sm","S+":"S-","Sm":"Sp","S-":"S+",
                "Sig":"SigDag","Tau":"TauDag","SigDag":"Sig","TauDag":"Tau"}

def get_dagger(self,conjugate=True):
    """Return the dagger of a multioperator"""
    # Same accumulation fix as jordan_wigner(): build the flat term list
    # directly instead of combining one term at a time with "+".
    outterms = []
    for opi in self.op: # loop over terms
        n = len(opi) # number of terms in the product
        c = opi[0] # coefficient
        newterm = [np.conjugate(c)] # initialize
        for i in range(n-1,0,-1): # reversed: (AB...)^dagger = ...B^dagger A^dagger
            name = opi[i][0] # name
            jj = opi[i][1] # index
            name2 = _dagger_name.get(name,name)
            newterm.append([name2,jj])
        outterms.append(newterm)
    out = MultiOperator(term=False)
    out.op = outterms
    out.i = len(outterms)-1
    return obj2MO(out) # return MO








