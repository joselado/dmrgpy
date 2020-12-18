import numpy as np


def imag2real(zs,us,x=None):
    """Translate from imaginary axis to real axis"""
    try: from .padetk import pade
    except: 
        print("Not functional yet")
        exit()
    if x is None: x = np.linspace(-4.,4.,300)
    p = pade.fit(zs, us)
    return x,np.array([p(ix) for ix in x])




