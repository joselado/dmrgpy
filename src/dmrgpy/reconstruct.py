import numpy as np
# routines for distribution reconstruction

try: from .maxenttk.pymaxent import reconstruct as reconstruct_from_moments
except:
    print("Reconstruction not functional yet")
    raise


def reconstruct_distribution(x,y,bnds=None,delta=1e-3,**kwargs):
    """Given a certain distribution, reconstruct it
    using the maximum entropy method"""
    yabs = np.abs(y)
    yabs = yabs/np.max(yabs) # maximum value
    if bnds is None: # autodetect bounds
        xnz = x[yabs>delta] # x with non-zero value
        bnds = [min(xnz),max(xnz)] # bounds of the distribution
    x0 = x[x<bnds[1]]
    y0 = y[x<bnds[1]]
    y0 = y0[x0>bnds[0]]
    x0 = x0[x0>bnds[0]]
    y0[y0<0.0] = 0.0 # set to zero
    print(bnds)
    xout,yout = reconstruct_distribution_restricted(x0,y0,**kwargs)
    from scipy.interpolate import interp1d
    y1 = interp1d(xout,yout,fill_value=0.0,bounds_error=False)(x)
    return x,y1

def reconstruct_distribution_restricted(x0,y0,n=6):
    """Given a certain distribution, reconstruct it
    using the maximum entropy method"""
    x,y = x0.real,y0.real # real part 
    # center the distribution in x=0
    y = y.real
    dx = x[1]-x[0] # interval
    scale0 = np.trapz(y,dx=dx) # scale for y axis
    y = y/scale0 # renormalize
    xmean = np.trapz(x*y,dx=dx) # average value
    x = x - xmean # shift the mean to zero
    # now compute the moments
    xscale = np.max(np.abs(x))
    xscale = np.sqrt(np.trapz(x**2*y,dx=dx)) # average value
#    print("x-scale",xscale)
    #xscale = np.max(np.abs(x[y/np.max(y)>1e-2])) # to interval 0,1
    x = x/xscale
    mu = [np.trapz(y*x**i,dx=dx) for i in range(n)] # compute all the moments
    scale = mu[0]
    mu = mu/scale # resscale
    bnds = [min(x)-dx,max(x)+dx]
    sol, lambdas = reconstruct_from_moments(mu,bnds=bnds)
    xo = x[x>bnds[0]]
    xo = xo[xo<bnds[1]]
    xout = xo*xscale+xmean # output x
    yout = sol(xo) # output y
    yout *= np.trapz(y0.real,dx=x0[1]-x0[0])/np.trapz(yout,dx=xout[1]-xout[0])
    return xout,yout


