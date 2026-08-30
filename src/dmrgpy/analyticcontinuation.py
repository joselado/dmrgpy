import numpy as np


def imag2real(zs,us,x=None):
    """Translate from imaginary axis to real axis"""
    try: 
        from .padetk import bruteforcepade as pade
#        from .padetk import taylorpade as pade
    except ImportError: 
        # this used to print "Not functional yet" and call exit(), i.e.
        # terminate the caller's whole process with status 0 -- from a
        # library, reachable straight from the public
        # get_dynamical_correlator(submode="CVMimag")
        raise NotImplementedError(
            "submode='CVMimag' needs the Pade analytic-continuation module "
            "dmrgpy.padetk, which is not part of this package; use "
            "submode='CVM'/'KPM'/'EX' instead")
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
    raise # not finished yet




