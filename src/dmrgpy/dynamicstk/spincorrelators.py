
import numpy as np

def get_full_SS_correlator(sc,mode="ED",i=0,j=None,**kwargs):
    """Return the full SS correlator, summing over ground states
    and summing over components. Factor out trivial contributions
    from the delta function"""
    if mode !="ED": raise # only implemented for ED
    if j is None: j = i # overwrite
    wfs = sc.get_gs_manifold(mode="ED") # return ground state manifold
    sc.get_gs(mode="ED") # return ground state manifold
    print("Computing SS correlator with ",len(wfs),"states")
    do = 0. # initialize
    for wf in wfs: # loop over wavefunctions
        sc.set_gs(wf) # set this wavefunction as ground state
        Si = [sc.Sx[i],sc.Sy[i],sc.Sz[i]] # spin operators
        Sj = [sc.Sx[j],sc.Sy[j],sc.Sz[j]] # spin operators
        for k in range(3): # loop over components
            Ai = Si[k]
            Aj = Sj[k]
#            Ai = Ai - wf.dot(Ai*wf) # remove background
#            Aj = Aj - wf.dot(Aj*wf) # remove background
            # compute correlator
            Iden = 4*Si[k]*Si[k] # identity
            (es,ds) = sc.get_dynamical_correlator(name=(Ai,Aj),
                    mode="ED",**kwargs)
            do += ds  # add contribution
    return es,do/(3*len(wfs)) # return result 





