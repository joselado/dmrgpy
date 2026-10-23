import numpy as np


def imag2real(zs,us,x=None):
    """Translate from imaginary axis to real axis"""
    from .padetk import bruteforcepade as pade
    p = pade.fit(zs,us)
    if x is None: x = np.linspace(-4.,4.,300)
    return x,np.array([p(ix) for ix in x])



def polyexp_fit(zs,us,n=10,bnds=[-5.,5.]):
    """Perform an analytic continuation using a polyexp fit"""
    def f(lamb):
      """Polynomial function"""
      return lambda x: np.exp(np.sum([lamb[i]*x**i for i in range(len(lamb))]))
    from scipy.integrate import quad
    def errorf(lamb):
        """Function to compute the error"""
        fx = f(lamb) # get the function to evaluate
        def gz(tau):
            """Do the transformation"""
            fint = lambda x: fx(x)*np.exp(-tau*x)/(1+np.exp(beta*x))
            out = quad(fint,bnds[0],bnds[1])
    from scipy.optimize import minimize
    raise NotImplementedError("polyexp_fit is not finished yet") # not finished yet




