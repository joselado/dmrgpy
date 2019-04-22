from .. import multioperator


def get_zero(name="multioperator"):
    return multioperator.MultiOperator(name,c=0.0) # generate the MO object

def get_si(j=0,**kwargs):
    if j==0: return get_sx(**kwargs)
    elif j==1: return get_sy(**kwargs)
    elif j==2: return get_sz(**kwargs)
    else: raise


def get_sx(name="multioperator",i=0):
    """Return the multioperator that computes sx in site i"""
    c = 1./2.
    mo = multioperator.MultiOperator(name,c=c) # generate the MO object
    mo.add_operator("Cdag",2*i+1) # add term
    mo.add_operator("C",2*i) # add term
    mo2 = multioperator.MultiOperator(name,c=c) # generate the MO object
#    mo.new_term(c=c) # generate a new term
    mo2.add_operator("Cdag",2*i) # add term
    mo2.add_operator("C",2*i+1) # add term
    return mo+mo2 # return multioperator



def get_sy(name="multioperator",i=0):
    """Return the multioperator that computes sx in site i"""
    c = 1j/2.
    mo = multioperator.MultiOperator(name,c=c) # generate the MO object
    mo.add_operator("Cdag",2*i+1) # add term
    mo.add_operator("C",2*i) # add term
    mo.new_term(c=-c) # generate a new term
    mo.add_operator("Cdag",2*i) # add term
    mo.add_operator("C",2*i+1) # add term
    return mo # return multioperator



def get_sz(name="multioperator",i=0):
    """Return the multioperator that computes sx in site i"""
    c = 1./2.
    mo = multioperator.MultiOperator(name,c=c) # generate the MO object
    mo.add_operator("Cdag",2*i) # add term
    mo.add_operator("C",2*i) # add term
    mo.new_term(c=-c) # generate a new term
    mo.add_operator("Cdag",2*i+1) # add term
    mo.add_operator("C",2*i+1) # add term
    return mo # return multioperator




def get_sysy(name="multioperator",i=0,j=0):
    return get_sy(name=name,i=i)*get_sy(name=name,i=j)


def get_szsz(name="multioperator",i=0,j=0):
    return get_sz(name=name,i=i)*get_sz(name=name,i=j)


def get_sxsx(name="multioperator",i=0,j=0):
    """Return the multioperator that computes SxSx in site i and j"""
    sxi = get_sx(name=name,i=i)
    sxj = get_sx(name=name,i=j)
    return sxi*sxj # return the product
    c = 1./4.
    # generate the MO object
    mo = multioperator.MultiOperator(name,c=c,term=False) 
    for ii in range(2): # loop over up/down
      for jj in range(2): # loop over up/down
        mo.new_term(c=c) # generate a new term
        mo.add_operator("Cdag",2*i+1-ii) # add term
        mo.add_operator("C",2*i+ii) # add term
        mo.add_operator("Cdag",2*j+1-jj) # add term
        mo.add_operator("C",2*j+jj) # add term
    return mo # return multioperator


def get_densitydensity_spinful(name="multioperator",i=0,j=0):
    """Return the multioperator that computes SxSx in site i and j"""
    c = 1./4.
    # generate the MO object
    mo = multioperator.MultiOperator(name,term=False) 
    for ii in range(2): # loop over up/down
      for jj in range(2): # loop over up/down
        mo.new_term(c=c) # generate a new term
        mo.add_operator("Cdag",2*i+ii) # add term
        mo.add_operator("C",2*i+ii) # add term
        mo.add_operator("Cdag",2*j+jj) # add term
        mo.add_operator("C",2*j+jj) # add term
    return mo # return multioperator






