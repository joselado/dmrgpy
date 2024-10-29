import numpy as np
from numba import jit
import scipy.linalg as lg
from numpy.polynomial.polynomial import polyfit
#from statsmodels.tsa.ar_model import sarimax 


import warnings

#import pmdarima as pm

warnings.filterwarnings("ignore")


def restrict(y):
    y2 = y+0.
    y2[y<-1.] = -1.
    y2[y>1.] = 1.
    return y2


def power_transform(mus):
    """Return the direct and inverse transformations for a power law"""
    def f(x): 
        scale = np.array(range(len(x))) +1.# scale
        return x*scale
    def finv(x): 
        scale = np.array(range(len(x))) +1.# scale
        return x/scale
    return f,finv # return functions

def no_transform():
    """Return the direct and inverse transformations for a power law"""
    return lambda x: x,lambda x: restrict(x)


def fit_power_transform(mus):
    """Fit to a power-law, and then transform"""
    y = np.abs(mus) # absolute value
    n = len(y) # number of points
    na = 3 # minimum number of points to average
    yi = [np.max(y[i:i+na]) for i in range(n-na)] # averge over several points
    yi = np.array(yi) # convert to array
    x = np.array(range(len(yi)))+1. # indexes
    yl = np.log(yi) # logarithm 
    xl = np.log(x) # logarithm
#    slope, intercept, r_value, p_value, std_err = linregress(xi,y)
    c,stats = polyfit(xl,yl,1,full=True,w=x) # fit to a power law
    def f(x): 
        scale = np.array(range(len(x))) +1.# scale
        scale = np.power(scale,c[1]) # scaling
        return x/np.exp(c[0])/scale
    def finv(x):
        x = x*np.exp(c[0]) # scale back
        scale = np.array(range(len(x))) +1.# scale
        scale = np.power(scale,c[1]) # scaling
        out = x*scale
        out = restrict(out)
        return out
    return f,finv





#def extrapolate_moments(mus,fac):
#    mus = extrapolate_moments_2(mus,np.sqrt(fac))
#    return extrapolate_moments_2(mus,np.sqrt(fac))

def extrapolate_moments(mus0,fac,extrapolation_mode="1/n"):
    """Extrapolate moments"""
    if np.max(mus0.imag)>1e-4: raise # not implemented
    mus = mus0.real
    if extrapolation_mode=="plain":
        ftrans,ftransinv = no_transform()
    elif extrapolation_mode=="1/n": 
        ftrans,ftransinv = power_transform(mus)
    elif extrapolation_mode=="power": 
        ftrans,ftransinv = fit_power_transform(mus)
    mus = ftrans(mus) # scale the moments
    L = len(mus)//2
    T = len(mus)
    L = T
    P = int(fac*T) # prediction
    train = mus[0:L].real # train data
    test = mus[L:T] # test data
#    model = AR(train).fit(ic="aic") # get the model
    lags = round(12*(len(train)/100.)**(1/4.))
    from statsmodels.tsa.ar_model import AR
    from statsmodels.tsa.arima_model import ARIMA
    from statsmodels.tsa.ar_model import AutoReg
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
    mus2 = np.zeros(P,dtype=np.complex128) 
    mus2[0:L] = mus[0:L] # initial data
    mus2[L:P] = pred[:] # predicted data
    mus2 = ftransinv(mus2) # transform back
#    print(extrapolation_mode,np.max(mus0),np.max(mus2))
    return mus2
#    mus2 = jackson_kernel(mus2) # use jackson kernel
#    return mus2






def deconvolution(es,gs,mode=None,delta0=1e-6,delta=None):
    if mode is None: return es,gs
    elif mode=="pm": # poor man deconvolution
        selfh = es - 1j/gs +1j*delta0 # finite size selfenergy
        self = es - 1j/gs - 1j*delta + 1j*delta0 # infinite selfenergy
        return es,1j/(es-self+1j*delta0) # return the sharpened one
    else: raise












