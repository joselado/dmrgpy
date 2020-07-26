import numpy as np
from numba import jit
import scipy.linalg as lg

from statsmodels.tsa.ar_model import AutoReg
#from statsmodels.tsa.ar_model import sarimax 


import warnings
from statsmodels.tsa.ar_model import AR

from statsmodels.tsa.arima_model import ARIMA

import pmdarima as pm

warnings.filterwarnings("ignore")


#def extrapolate_moments(mus,fac):
#    mus = extrapolate_moments_2(mus,np.sqrt(fac))
#    return extrapolate_moments_2(mus,np.sqrt(fac))

def extrapolate_moments(mus,fac):
    """Extrapolate moments"""
    L = len(mus)//2
    T = len(mus)
    L = T
    P = int(fac*T) # prediction
    train = mus[0:L].real # train data
    test = mus[L:T] # test data
#    model = AR(train).fit(ic="aic") # get the model
    lags = round(12*(len(train)/100.)**(1/4.))
    model = AutoReg(train,lags=lags,trend="ct").fit(cov_type="HC1") # get the model

#    model = pm.auto_arima(train, start_p=1, start_q=1,
#                         test='adf',
#                         max_p=3, max_q=3, m=10,
#                         start_P=0, seasonal=True,
#                         d=None, D=1, trace=True,
#                         error_action='ignore',  
#                         suppress_warnings=True, 
#                         stepwise=True)


#    pred = model.predict(n_periods=P-L) # prediction
    pred = model.predict(start=L,end=P-1) # prediction
    mus2 = np.zeros(P,dtype=np.complex) 
    mus2[0:L] = mus[0:L] # initial data
    mus2[L:P] = pred[:] # predicted data
    return mus2
#    mus2 = jackson_kernel(mus2) # use jackson kernel
#    return mus2

