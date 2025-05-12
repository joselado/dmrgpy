# function to deal with fermionic operators in very long chains

from ..multioperator import MultiOperator,jordan_wigner

def is_long_operator(MO):
    """Check if we are dealing with a long operator"""
    pass


def toMPO(MBO,MO,maxp=4):
    """Transform a long operator into an MPO by using products"""
    if not type(MO)==MultiOperator: raise
    if len(MO.op)>1: raise # not implemented
    MO = jordan_wigner(MO) # use JW if needed
    opt = MO.op[0] # total operator
    nump = len(opt) - 1 # total number of products
    opsi = [] # empty list of the new products
    for i in range(nump): # loop over products
        if i%maxp==0: # create a new one
            if i==0: 
                opi = [opt[0]] # store coefficient
            else: 
                opsi.append(opi) # store previous one
                opi = [1.0] # default coefficient
        opi.append(opt[i+1]) # store
    if len(opi)>1: # store last one
        opsi.append(opi) # store last one
    MO0 = MO.copy() # make a copy
    out = None # start
#    print(opt)
    for opi in opsi: # loop over terms
#        print(opi)
        MOi = MO0.copy() # make a copy
        MOi.op = [opi] # set the term
        if out is None:
            out = MBO.toMPO(MOi) # add contribution
        else:
            out = out*MBO.toMPO(MOi) # add contribution
    return out # return result



