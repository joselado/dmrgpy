import numpy as np
from scipy.optimize import minimize

def fit(f,x0,ntries=10,method="Powell"):
    """Fit a function several times, and estimate the error"""
    error = 1e8
    x = 0
    for i in range(ntries):
        opt = minimize(f,x0,method=method) # minimize function
        if f(opt.x)<error: 
            x = opt.x + 0. # store
            error = f(opt.x)
        x0 = np.random.random(x0.shape) # new guess
    return x



