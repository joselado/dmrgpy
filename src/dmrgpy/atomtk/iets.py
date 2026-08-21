import numpy as np


def check_iorb(iorb,norb):
    """Check that a requested orbital index is in range.

    Both entry points below used to overwrite their caller's `iorb` with
    a hardcoded 0 immediately after unpacking it, so every call silently
    computed orbital 0 no matter what was asked for -- a factor ~2 in
    spectral weight and a completely different peak position for a d-shell
    atom, with no warning. Those two lines are gone; this check exists so
    that an out-of-range `iorb` raises instead of turning into a fresh
    variant of the same silently-wrong behaviour."""
    if not (0<=iorb<norb):
        raise ValueError("iorb must be in [0,%d), got %d"%(norb,iorb))


def get_orbital_cotunneling(fc,es=[0.],mode="ED",
                            submode="ED",delta=1e-3,
                            iorb=0,
                            **kwargs):
    """Get the orbital cotunneling component of the IETS"""
    es = np.array(es) ; es = es[es>0.] # only positive energies
    ne = len(es) # number of energies
    
    wf = fc.get_gs(mode=mode)

    # create fast operators
    Cdagup = [fc.toMPO(A,mode=mode) for A in fc.Cdagup] # create fast operator
    Cup = [fc.toMPO(A,mode=mode) for A in fc.Cup] # create fast operator
    Cdagdn = [fc.toMPO(A,mode=mode) for A in fc.Cdagdn] # create fast operator
    Cdn = [fc.toMPO(A,mode=mode) for A in fc.Cdn] # create fast operator

    y_neg, y_pos = 0.,0. # initialize

    norb = len(Cup) # number of orbitals actually present in this chain
    check_iorb(iorb,norb) # sanity check on the requested input orbital
    for iorbj in range(norb): # loop over output orbitals
        if iorb==iorbj: continue # skip this iteration
        for Op1 in [Cdagup[iorb],Cdagdn[iorb]]: # loop over input orbital
            for Op2 in [Cup[iorbj],Cdn[iorbj]]: # loop over output orbital
                A01 = Op1*Op2 # in and out
                A01 = A01 - wf.dot(A01*wf) # minus the average 
                T01 = (A01,A01.get_dagger()) # positive bias operator
                T10 = (A01.get_dagger(),A01) # negative bias operator
                x_pos,ypt = fc.get_dynamical_correlator(name=T01,
                                                        mode=mode,submode=submode,es=es,delta=delta,**kwargs)
                x_neg,ynt = fc.get_dynamical_correlator(name=T10,
                                                        mode=mode,submode=submode,es=es,delta=delta,**kwargs)
                y_pos = y_pos + ypt # add contribution
                y_neg = y_neg + ynt # add contribution
        
    x_neg = -x_neg # invert values
    
    # Sort values for negative energies
    sorted_indices_neg = np.argsort(x_neg)
    x_neg = x_neg[sorted_indices_neg]
    y_neg = y_neg[sorted_indices_neg]
    
    # Combine positive and negative
    x_combined = np.concatenate([x_neg, x_pos]) 
    y_combined = np.concatenate([y_neg, y_pos])
    
    return x_combined, y_combined    


def get_spinflip(fc,es=[0.],mode="ED",submode="ED",delta=1e-3,
                 iorb=0, # orbital for the spin flip
                 **kwargs,
                 ):
    """Get the spin flip excitations"""
    es = np.array(es) ; es = es[es>0.] # only positive energies
    ne=len(es) # number of energies
 
    check_iorb(iorb,len(fc.Sx)) # sanity check on the requested orbital
    Sis = [fc.Sx[iorb],fc.Sy[iorb],fc.Sz[iorb]] # spin operators in the input orbital
    y_sum = 0. # initialize
    for Si in Sis: # loop over spin operators
        name = (Si - fc.vev(Si,mode=mode), Si - fc.vev(Si,mode=mode)) # correlator to compute
        x, y = fc.get_dynamical_correlator(delta=delta, name=name, 
                es=es,mode=mode,submode=submode,**kwargs)
        y_sum = y_sum + y # add contribution
    
    # add the negative energies
    y_neg = y_sum
    x_neg = -x
    
    sorted_indices_neg = np.argsort(x_neg)
    x_neg = x_neg[sorted_indices_neg]
    y_neg = y_neg[sorted_indices_neg]
    y_accum = np.concatenate([y_neg, y_sum])
    
    x_array = np.concatenate([x_neg,x])
    
    return x_array, y_accum
