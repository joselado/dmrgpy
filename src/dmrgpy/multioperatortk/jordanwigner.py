from .. import multioperator

def obj2MO(a): return multioperator.obj2MO([a])

def CdagC(i,j):
    if i==j: return obj2MO(["Adag",i])*obj2MO(["A",j])
    elif i<j:
        m = obj2MO(["Adag",i])
        for k in range(i,j-1):
            m = m*obj2MO(["F",k+1])
        return m*obj2MO(["A",j])
    elif j<i: return -1*CCdag(j,i)


def CCdag(i,j):
    if i==j: return obj2MO(["A",i])*obj2MO(["Adag",j])
    elif i<j:
        m = -1*obj2MO(["A",i])
        for k in range(i,j-1):
            m = m*obj2MO(["F",k+1])
        return m*obj2MO(["Adag",j])
    elif j<i: return -1*CdagC(j,i)


def C(i):
    m = 1.0
    for k in range(i):
        m = obj2MO(["F",k])*m
    return m*obj2MO(["A",i])


def Cdag(i):
    m = obj2MO(["Adag",i])
    for k in range(i):
        m = obj2MO(["F",k])*m
    return m


def CC(i,j):
    if i==j: return 0
    elif i<j:
        m = -1*obj2MO(["A",i])
        for k in range(i,j-1):
            m = m*obj2MO(["F",k+1])
        return m*obj2MO(["A",j])
    elif j<i: return -1*CC(j,i)


def CdagCdag(i,j):
    if i==j: return 0
    elif i<j:
        m = -1*obj2MO(["Adag",i])
        for k in range(i,j-1):
            m = m*obj2MO(["F",k+1])
        return m*obj2MO(["Adag",j])
    elif j<i: return -1*CdagCdag(j,i)

def CdagCCdagC(i,j,k,l):
    return CdagC(i,j)*CdagC(k,l)

def one_fermion(ni,i,**kwargs):
    if (ni=="C"): return C(i,**kwargs)
    elif (ni=="Cdag"): return Cdag(i,**kwargs)
    else: raise


def two_fermions(ni,i,nj,j,**kwargs):
    if (ni=="C") and (nj=="Cdag"): return CCdag(i,j,**kwargs)
    elif (ni=="Cdag") and (nj=="C"): return CdagC(i,j,**kwargs)
    elif (ni=="Cdag") and (nj=="Cdag"): return CdagCdag(i,j,**kwargs)
    elif (ni=="C") and (nj=="C"): return CC(i,j,**kwargs)
    else: raise


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



