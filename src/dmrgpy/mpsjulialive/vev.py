

def vev(MBO,MO):
    """Compute a vacuum expectation value"""
    wf0 = MBO.get_gs() # get ground state
    return wf0.dot(MO*wf0) # return the expectation value


