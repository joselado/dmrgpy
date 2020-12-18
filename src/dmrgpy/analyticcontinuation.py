import numpy as np


def imag2real(zs,us,x=None):
    """Translate from imaginary axis to real axis"""
    try: 
        from .padetk import bruteforcepade as pade
#        from .padetk import taylorpade as pade
    except: 
        print("Not functional yet")
        exit()
    p = pade.fit(zs,us)
    if x is None: x = np.linspace(-4.,4.,300)
    return x,np.array([p(ix) for ix in x])




