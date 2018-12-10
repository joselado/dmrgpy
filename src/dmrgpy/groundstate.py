


def best_gs(sc,n=1):
    """Compute many ground states, and retain only the best one"""
    wfs = [] # list with GS
    es = [] # list with energies
    emin = 1e8 # grund state energy
    wf0 = None
    for i in range(n): # loop
        sc.computed_gs = False # initialize
        e0 = sc.gs_energy() # gounrd state energy
        print(e0)
        if e0<emin: wf0 = sc.wf0 # copy wavefunction
    sc.set_initial_wf(wf0)


