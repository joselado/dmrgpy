from .. import multioperator

def get_sx(name="multioperator",i=0):
    """Return the multioperator that computes sx in site i"""
    c = 1./2.
    mo = multioperator.MultiOperator(name,c=c) # generate the MO object
    mo.add_operator("Cdag",2*i+1) # add term
    mo.add_operator("C",2*i) # add term
    mo.new_term(c=c) # generate a new term
    mo.add_operator("Cdag",2*i) # add term
    mo.add_operator("C",2*i+1) # add term
    return mo # return multioperator



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







def get_sxsx(name="multioperator",i=0,j=0):
    """Return the multioperator that computes SxSx in site i and j"""
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






