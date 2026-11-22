from .. import multioperator

def obj2MO(a): return multioperator.obj2MO([a])


def _term(c, factors):
    """The one-term MultiOperator c*factors[0]*factors[1]*..., built in
    one step. Building a string by one product per F factor, as this
    module used to, copies the growing term at every step, so a term of
    span s cost O(s^2); a rotated hopping of a Gaussian basis spans the
    chain, and its expansion dominated the build of the Hamiltonian. The
    terms are the same as those products give, coefficient and factor
    order included (tests/test_jordan_wigner_strings.py)."""
    mo = multioperator.MultiOperator(name="multioperator", term=False)
    mo.op = [[c]+[[name, k] for (name, k) in factors]]
    mo.i = 0
    return mo


def _string(i, j):
    """The F factors strictly between sites i < j."""
    return [("F", k) for k in range(i+1, j)]


def CdagC(i,j):
    if i==j: return _term(1.0, [("Adag",i), ("A",j)])
    elif i<j: return _term(1.0, [("Adag",i)]+_string(i,j)+[("A",j)])
    else: return _term(1.0, [("A",j)]+_string(j,i)+[("Adag",i)])


def CCdag(i,j):
    if i==j: return _term(1.0, [("A",i), ("Adag",j)])
    elif i<j: return _term(-1.0, [("A",i)]+_string(i,j)+[("Adag",j)])
    else: return _term(-1.0, [("Adag",j)]+_string(j,i)+[("A",i)])


def C(i):
    return _term(1.0, [("F",k) for k in range(i-1,-1,-1)]+[("A",i)])


def Cdag(i):
    return _term(1.0, [("F",k) for k in range(i-1,-1,-1)]+[("Adag",i)])


def CC(i,j):
    # C_i C_i = 0 by Pauli exclusion, but it has to be returned as an
    # identically-zero MultiOperator rather than the Python int 0: every
    # caller multiplies this into a term and then asks for .op, so the int
    # surfaced as "AttributeError: 'float' object has no attribute 'op'"
    # on all three DMRG backends for a perfectly legitimate vev(C[i]*C[i])
    # (mode="ED" answered 0j). The product below is the same form the
    # Cdag path already produces, and every backend evaluates it to 0.
    if i==j: return _term(1.0, [("A",i), ("A",i)])
    elif i<j: return _term(-1.0, [("A",i)]+_string(i,j)+[("A",j)])
    else: return _term(1.0, [("A",j)]+_string(j,i)+[("A",i)])


def CdagCdag(i,j):
    return CC(j,i).get_dagger()

def CdagCCdagC(i,j,k,l):
    return CdagC(i,j)*CdagC(k,l)

def one_fermion(ni,i,**kwargs):
    if (ni=="C"): return C(i,**kwargs)
    elif (ni=="Cdag"): return Cdag(i,**kwargs)
    else: raise ValueError("Unknown fermion operator name %s" % ni)


def two_fermions(ni,i,nj,j,**kwargs):
    if (ni=="C") and (nj=="Cdag"): return CCdag(i,j,**kwargs)
    elif (ni=="Cdag") and (nj=="C"): return CdagC(i,j,**kwargs)
    elif (ni=="Cdag") and (nj=="Cdag"): return CdagCdag(i,j,**kwargs)
    elif (ni=="C") and (nj=="C"): return CC(i,j,**kwargs)
    else: raise ValueError("Unknown fermion operator names %s, %s" % (ni,nj))


def four_fermions(ni,i,nj,j,nk,k,nl,l,**kwargs):
    return two_fermions(ni,i,nj,j,**kwargs)*two_fermions(nk,k,nl,l,**kwargs)


def product2jw(ns,inds):
    """Transform a generic product of operators into bosonic
    using Jordan Wigner, the transformation may not be optimal"""
    out = 1
    for (ni,i) in zip(ns,inds):
        if ni=="C" or ni=="Cdag": out = out*one_fermion(ni,i)
        else: out = out*obj2MO([ni,i])
    return out



