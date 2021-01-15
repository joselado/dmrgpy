import numpy as np
import scipy.linalg as lg
from scipy.stats import linregress
from numpy.polynomial.polynomial import polyfit


def size_extrapolation(x,y,n=1):
    """Perform a size extrapolation using a power law
    Y(X) = A + B/X"""
    xi = (1./np.array(x)) # inverse
#    slope, intercept, r_value, p_value, std_err = linregress(xi,y)
    c,stats = polyfit(xi,y,n,full=True,w=np.array(x))
    return c[0]



#from sklearn.decomposition import PCA
#pca = PCA(n_components=2)
#pca.fit(X)


