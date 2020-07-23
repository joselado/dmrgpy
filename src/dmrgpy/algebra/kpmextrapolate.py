import numpy as np
from numba import jit
import scipy.linalg as lg

from statsmodels.tsa.ar_model import AutoReg
#from statsmodels.tsa.ar_model import sarimax 


import warnings
from statsmodels.tsa.ar_model import AR

warnings.filterwarnings("ignore")

def extrapolate_moments(mus,fac):
    """Extrapolate moments"""
    L = len(mus)//2
    T = len(mus)
    L = T
    P = int(fac*T) # prediction
    train = mus[0:L].real # train data
    test = mus[L:T] # test data
    model = AR(train).fit(ic="aic") # get the model
#    model = AutoReg(train,lags=2,trend="c").fit() # get the model
#    model = SARIMAX(train).fit() # get the model
#    model = sarimax.SARIMAX(train).fit() # get the model
    pred = model.predict(start=L,end=P-1) # prediction
#    error = test - model.predict(start=1,end=L-1)
#    print("Error",np.mean(np.abs(error)))
    mus2 = np.zeros(P,dtype=np.complex) 
    mus2[0:L] = mus[0:L] # initial data
    mus2[L:P] = pred[:] # predicted data
    return mus2
#    mus2 = jackson_kernel(mus2) # use jackson kernel
#    return mus2

