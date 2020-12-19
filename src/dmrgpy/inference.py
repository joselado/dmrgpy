import numpy as np

def points2function(x,y,bounds=None):
    """Given a set of points for a holomorphic function,
    derive the continuum function that fullfills the
    analiticity condition"""
    yr = y.real # real part
    yi = y.imag # imaginary part
    from scipy.signal import hilbert
    from scipy.interpolate import interp1d
    if bounds is None: bounds = [min(x),max(x)]
    xc = np.linspace(bounds[0],bounds[1],80)
#    xc = x
    def f(yc):
        """Function to minimize"""
        yrc = interp1d(xc,yc,fill_value=0.0,bounds_error=False)(x) # predicted
#        yic = interp1d(xc,hilbert(yc).imag,fill_value=0.0,bounds_error=False)(x) # predicted
        out = np.sum((yrc-yr)**2) #+ np.sum((yic-yi)**2)
        print(out)
        return out
    from scipy.optimize import minimize
    bndx = [(0,10) for ix in xc]
    bndx = None
    error = 1e3
    ntries = 30
    for i in range(ntries):
      y0 = np.random.random(xc.shape[0])
      res = minimize(f,y0,method="SLSQP",bounds=bndx,
            options={"ftol":1e-10,"maxiter":1000})
      if f(res.x)<error:
          outy = res.x
          error = f(res.x)
    return xc,hilbert(outy) # return optimized distribution


