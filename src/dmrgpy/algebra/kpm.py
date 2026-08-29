# kernel polynomial method libraries
from __future__ import print_function,division
import scipy.sparse.linalg as lg
from scipy.sparse import csc_matrix as csc
import numpy.random as rand
from scipy.sparse import coo_matrix,csc_matrix,bmat
import numpy as np
from scipy.signal import hilbert
from . import algebra

# check that the fortran library exists
try: 
  from . import kpmf90 
  use_fortran = True
except:
  use_fortran = False # use python routines
#  print("FORTRAN library not present, using default python one")




def get_moments(v,m,n=100,use_fortran=use_fortran,test=False):
  """ Get the first n moments of a certain vector
  using the Chebychev recursion relations"""
  if use_fortran:
    from .kpmf90 import get_momentsf90 # fortran routine
    mo = coo_matrix(m) # convert to coo matrix
    v = algebra.matrix2vector(v)
# call the fortran routine
    mus = get_momentsf90(mo.row+1,mo.col+1,mo.data,v,n) 
    return mus # return fortran result
  else:
   if test: return python_kpm_moments_clear(v,m,n=n)
   else: return python_kpm_moments(v,m,n=n)


def python_kpm_moments(v,m,n=100):
  """Python routine to calculate moments"""
  mus = np.array([0.0j for i in range(2*n)]) # empty arrray for the moments
  am = v.copy() # zero vector
  a = m@v  # vector number 1
  bk = algebra.braket_ww(v,v)
#  bk = (np.transpose(np.conjugate(v))*v)[0,0] # scalar product
  bk1 = algebra.braket_ww(a,v)
#  bk1 = (np.transpose(np.conjugate(a))*v)[0,0] # scalar product
  
  mus[0] = bk.copy()  # mu0
  mus[1] = bk1.copy() # mu1
  for i in range(1,n): 
    ap = 2*m@a - am # recursion relation
    bk = algebra.braket_ww(a,a)
#    bk = (np.transpose(np.conjugate(a))*a)[0,0] # scalar product
    bk1 = algebra.braket_ww(ap,a)
#    bk1 = (np.transpose(np.conjugate(ap))*a)[0,0] # scalar product
    mus[2*i] = 2.*bk
    mus[2*i+1] = 2.*bk1
    am = a.copy() # new variables
    a = ap.copy() # new variables
  mu0 = mus[0] # first
  mu1 = mus[1] # second
  for i in range(1,n): 
    mus[2*i] +=  - mu0
    mus[2*i+1] += -mu1 
  return mus


def python_kpm_moments_clear(v,m,n=100):
  """Python routine to calculate moments"""
  mus = np.array([0.0j for i in range(2*n)]) # empty arrray for the moments
  a0 = v.copy() # first vector
  am = v.copy() # first vector
  a = m*v  # vector number 1
  mus[0] = 1.  # mu0
  mu = (np.transpose(np.conjugate(a0))*a)[0,0] # scalar product
  mus[1] = mu # mu1
  for i in range(1,2*n): 
    ap = 2*m*a - am # recursion relation
    mu = (np.transpose(np.conjugate(a0))*a)[0,0] # scalar product
    mus[i] = mu # store
    am = a.copy() # new variables
    a = ap.copy() # new variables
  return mus






def get_momentsA(v,m,n=100,A=None):
  """ Get the first n moments of a certain vector
  using the Chebychev recursion relations"""
  mus = np.array([0.0j for i in range(n)]) # empty arrray for the moments
  am = algebra.matrix2vector(v) # zero vector
  a = m@v  # vector number 1
#  print(v.shape,A.shape)
  bk = algebra.braket_wAw(v,A,v)
#  bk = (np.transpose(np.conjugate(v))*A*v)[0,0] # scalar product
  bk1 = algebra.braket_wAw(a,A,v)
#  bk1 = (np.transpose(np.conjugate(a))*A*v)[0,0] # scalar product
  mus[0] = bk  # mu0
  mus[1] = bk1 # mu1
  for i in range(2,n): 
#    print(A)
    ap = 2.*m@a - am # recursion relation
    bk = algebra.braket_wAw(ap,A,v)
#    bk = (np.transpose(np.conjugate(ap))*A*v)[0,0] # scalar product
    mus[i] = bk
    am = a.copy() # new variables
    a = ap.copy() # new variables
  mu0 = mus[0] # first
  mu1 = mus[1] # second
  return mus




def get_moments_ij(m0,n=100,i=0,j=0,use_fortran=use_fortran):
  """ Get the first n moments of a the |i><j| operator
  using the Chebychev recursion relations"""
  m = coo_matrix(m0,dtype=np.complex128)
  if use_fortran:
    mus = kpmf90.get_moments_ij(m.row+1,m.col+1,m.data,n,m.shape[0],i+1,j+1)
    return mus
  else:
    mus = np.zeros(n,dtype=np.complex128) # empty arrray for the moments
    v = np.zeros(m.shape[0],dtype=np.complex128) ; v[i] = 1.0 # initial vector
    v = np.matrix([v]).T # zero vector
    am = v.copy()
    a = m*v  # vector number 1
    bk = v[j] # scalar product
    bk1 = a[j,0] # scalar product
    mus[0] = bk  # mu0
    mus[1] = bk1 # mu1
    for ii in range(2,n): 
      ap = 2.*m*a - am # recursion relation
      bk = ap[j,0] # scalar product
      mus[ii] = bk
      am = a.copy() # new variables
      a = ap.copy() # new variables
    return mus

def get_moments_vivj(m0,vi,vj,n=100,use_fortran=False):
  if not use_fortran: return get_moments_vivj_python(m0,vi,vj,n=n)
  else: return get_moments_vivj_fortran(m0,vi,vj,n=n)


def get_moments_vivj_python(m0,vi,vj,n=100):
  """ Get the first n moments of a the |i><j| operator
  using the Chebychev recursion relations"""
  m = csc_matrix(m0,dtype=np.complex128)
  mus = np.zeros(n,dtype=np.complex128) # empty arrray for the moments
  v = vi.copy()
  am = v.copy()
  a = m@v  # vector number 1
  bk = algebra.braket_ww(vj,v)
#  bk = (vj.H*v).todense().trace()[0,0] # calculate bk
  bk1 = algebra.braket_ww(vj,a)
#  bk1 = (vj.H*a).todense().trace()[0,0] # calculate bk
  mus[0] = bk  # mu0
  mus[1] = bk1 # mu1
  for ii in range(2,n): 
    ap = 2.*m@a - am # recursion relation
    bk = algebra.braket_ww(vj,ap)
#    bk = (vj.H*ap).todense().trace()[0,0]
    mus[ii] = bk
    am = a.copy() # new variables
    a = ap.copy() # new variables
  return mus



def get_moments_vivj_fortran(m0,vi,vj,n=100):
    raise # I haven't check this function
    mo = coo_matrix(m0) # convert to coo matrix
    vi1 = vi.todense() # convert to conventional vector
    vj1 = vj.todense() # convert to conventional vector
# call the fortran routine
    mus = get_moments_vivj(mo.row+1,mo.col+1,mo.data,vi,vj,n)
    return mus # return fortran result



def get_mu_n_nh(hs,hs_dag,wfa,wfb,n=100):
  """Non-Hermitian KPM biorthogonal Chebyshev moments (port of NHKPM.jl's
  get_vn_NH/get_mu_n_NH, https://github.com/GUANGZECHEN/NHKPM.jl).

  hs, hs_dag are the already shifted-and-rescaled operator (z*I-H)/E_max
  and its conjugate transpose (matrices, not necessarily Hermitian - hs_dag
  need not equal hs). wfa,wfb are the ket/bra-side wavefunctions (already
  dressed by whichever operators define the correlator). Returns mu_n,
  length n, complex, with mu_n[k] = <wfb|vn[2k-1]> in the notation of the
  reference (only the odd-order Chebyshev vectors of the coupled H/H^dagger
  recursion enter the final reconstruction, see spec_from_moments_nh)."""
  hs = csc_matrix(hs,dtype=np.complex128)
  hs_dag = csc_matrix(hs_dag,dtype=np.complex128)
  v = algebra.matrix2vector(wfa)
  wfb = algebra.matrix2vector(wfb)

  alpha_prev2 = hs_dag@v                                  # alpha[1]
  alpha_prev1 = 2.*(hs@alpha_prev2) - v                   # alpha[2]

  vn_prev2 = v                                            # vn[1]
  vn_prev1 = 2.*(hs@vn_prev2)                             # vn[2]

  mus = np.zeros(n,dtype=np.complex128)
  mus[0] = algebra.braket_ww(wfb,vn_prev2)
  for k in range(1,n):
    alpha_x = 2.*(hs_dag@alpha_prev1) - alpha_prev2
    alpha_y = 2.*(hs@alpha_x) - alpha_prev1
    vn_x = 2.*alpha_prev1 + 2.*(hs_dag@vn_prev1) - vn_prev2
    vn_y = 2.*(hs@vn_x) - vn_prev1

    mus[k] = algebra.braket_ww(wfb,vn_x)

    alpha_prev2,alpha_prev1 = alpha_x,alpha_y
    vn_prev2,vn_prev1 = vn_x,vn_y
  return mus


def spec_from_moments_nh(mu_n,kernel="jackson"):
  """Reconstruct the non-Hermitian KPM spectral function rho(z) from the
  moments mu_n returned by get_mu_n_nh (see get_spec_kpm_nh), following
  NHKPM.jl's get_spec_kpm_NH: rho = sum_k gn[2k-1]*mu_n[k]*sin((2k-1)*pi/2),
  with gn the requested kernel's coefficients of order 2*len(mu_n)-1."""
  n = len(mu_n)
  ones = np.ones(2*n-1,dtype=np.complex128)
  if kernel=="jackson": gn_full = jackson_kernel(ones)
  elif kernel=="lorentz": gn_full = lorentz_kernel(ones)
  elif kernel in ("dirichlet","plain",None): gn_full = ones
  else: raise ValueError("Unknown kernel: "+str(kernel))
  gn = gn_full[0::2] # odd Chebyshev orders only (1-based indices 1,3,5,...)
  signs = (-1.)**np.arange(n) # sin((2k-1)*pi/2) for 0-based k
  return np.sum(gn*np.asarray(mu_n)*signs)


def get_spec_kpm_nh(h,z,wfa,wfb,E_max,n=100,kernel="jackson"):
  """Non-Hermitian KPM dynamical correlator/spectral function at a single
  complex frequency z: <wfb| delta-like-kernel(z-h) |wfa>, computed by
  rescaling hs=(z*I-h)/E_max and running the biorthogonal Chebyshev
  recursion (get_mu_n_nh) followed by the kernel reconstruction
  (spec_from_moments_nh). Unlike ordinary (Hermitian) KPM, the moments
  depend on z itself, so this recomputes them from scratch for every z -
  matching the reference algorithm's own cost profile."""
  dim = h.shape[0]
  hs = (-csc_matrix(h,dtype=np.complex128)+z*csc_matrix(
          (np.ones(dim),(range(dim),range(dim))),shape=(dim,dim),
          dtype=np.complex128))/E_max
  hs_dag = algebra.dagger(hs.todense())
  mu_n = get_mu_n_nh(hs,hs_dag,wfa,wfb,n=n)
  return spec_from_moments_nh(mu_n,kernel=kernel)



def full_trace(m_in,n=200,use_fortran=use_fortran):
  """ Get full trace of the matrix"""
  m = csc(m_in) # saprse matrix
  nd = m.shape[0] # length of the matrix
  mus = np.array([0.0j for i in range(2*n)])
#  for i in range(ntries):
  for i in range(nd):
    mus += local_dos(m_in,i=i,n=n,use_fortran=use_fortran)
  return mus/nd









def local_dos(m_in,i=0,n=200,use_fortran=use_fortran):
  """ Calculates local DOS using the KPM"""
  m = csc(m_in) # saprse matrix
  nd = m.shape[0] # length of the matrix
  mus = np.array([0.0j for j in range(2*n)])
  v = rand.random(nd)*0.
  v[i] = 1.0 # vector only in site i 
  v = csc(v).transpose()
# get the chebychev moments
  mus += get_moments(v,m,n=n,use_fortran=use_fortran) 
  return mus



def ldos(m_in,i=0,scale=10.,npol=None,ne=500,kernel="jackson"):
  """Return two arrays with energies and local DOS"""
  if npol is None: npol = ne
  mus = local_dos(csc_matrix(m_in)/scale,i=i,n=npol) # get coefficients
  xs = np.linspace(-1.0,1.0,ne,endpoint=True)*0.99 # energies
  ys = generate_profile(mus,xs,kernel=kernel)
  return (scale*xs,ys/scale)



ldos0d = ldos



def tdos(m_in,scale=10.,npol=None,ne=500,kernel="jackson",
              ntries=20,ewindow=None,frand=None,
              operator=None):
  """Return two arrays with energies and local DOS"""
  if npol is None: npol = ne
  mus = random_trace(m_in/scale,ntries=ntries,n=npol,fun=frand,
          operator=operator) 
  if ewindow is None or abs(ewindow)>scale: # no window provided
    xs = np.linspace(-1.0,1.0,ne,endpoint=True)*0.99 # energies
  else:
    xx = abs(ewindow/scale) # scale
    xs = np.linspace(-xx,xx,ne,endpoint=True)*0.99 # energies
  ys = generate_profile(mus,xs,kernel=kernel).real
  return (scale*xs,ys/scale)


tdos0d = tdos # redefine


def total_energy(m_in,scale=10.,npol=None,ne=500,ntries=20):
   x,y = tdos0d(m_in,scale=scale,npol=npol,ne=ne,ntries=ntries)
   z = .5*(np.sign(x)+1.)*x*y # function to integrate
   return np.trapezoid(z,x)



def random_trace(m_in,ntries=20,n=200,fun=None,operator=None):
  """ Calculates local DOS using the KPM"""
  if fun is not None: # check that dimensions are fine
    v0 = fun()
    if len(v0) != m_in.shape[0]: raise
  if fun is None:
#    def fun(): return rand.random(nd) -.5 + 1j*rand.random(nd) -.5j
    def fun(): return (rand.random(nd) - 0.5)*np.exp(2*1j*np.pi*rand.random(nd))
  m = csc(m_in) # saprse matrix
  nd = m.shape[0] # length of the matrix
  def pfun(x):
    v = fun()
    v = v/np.sqrt(v.dot(np.conjugate(v))) # normalize the vector
#    v = csc(v).transpose()
    if operator is None:
      mus = get_moments(v,m,n=n) # get the chebychev moments
    else:
      mus = get_momentsA(v,m,n=2*n,A=operator) # get the chebychev moments
    return mus
#  from . import parallel
#  out = [pfun(i) for i in range(ntries)] # perform all the computations
  from . import parallel
  out = parallel.pcall(pfun,range(ntries))
  mus = np.zeros(out[0].shape,dtype=np.complex128)
  for o in out: mus = mus + o # add contribution
  return mus/ntries



def random_trace_A(m_in,ntries=20,n=200,A=None):
  """ Calculates local DOS using the KPM"""
  m = csc(m_in) # saprse matrix
  nd = m.shape[0] # length of the matrix
  mus = np.array([0.0j for j in range(n)])
  for i in range(ntries): # loop over tries
    #v = rand.random(nd) - .5
    v = rand.random(nd) -.5 + 1j*rand.random(nd) -.5j
    v = v/np.sqrt(v.dot(v)) # normalize the vector
    v = csc(v).transpose()
    mus += get_momentsA(v,m,n=n,A=A) # get the chebychev moments
  return mus/ntries



def full_trace_A(m_in,ntries=20,n=200,A=None):
  """ Calculates local DOS using the KPM"""
  m = csc(m_in) # saprse matrix
  nd = m.shape[0] # length of the matrix
  mus = np.array([0.0j for j in range(2*n)])
  for i in range(nd): # loop over tries
    #v = rand.random(nd) - .5
    v = rand.random(nd)*0.
    v[i] = 1.0 # vector only in site i 
    v = csc(v).transpose()
    mus += get_momentsA(v,m,n=n,A=A) # get the chebychev moments
  return mus/nd



def correlator0d(m_in,i=0,j=0,scale=10.,npol=None,ne=500,write=True,
  x=None):
  """Return two arrays with energies and local DOS"""
  if npol is None: npol = ne
  mus = get_moments_ij(m_in/scale,n=npol,i=i,j=j,use_fortran=True)
  if np.sum(np.abs(mus.imag))>0.001:
#    print("WARNING, off diagonal has nonzero imaginary elements",np.sum(np.abs(mus.imag)))
     pass
  if x is None: xs = np.linspace(-1.0,1.0,ne,endpoint=True)*0.99 # energies
  else: xs = x/scale # use from input
  ys = generate_green_profile(mus,xs,kernel="jackson",use_fortran=False)/scale*np.pi # so it is the Green function
#  imys = hilbert(ys).imag
  if write: np.savetxt("CORRELATOR_KPM.OUT",np.matrix([scale*xs,-ys.imag,ys.real]).T)
  return (scale*xs,ys.real,ys.imag)




def dm_ij_energy(m_in,i=0,j=0,scale=10.,npol=None,ne=500,x=None):
  """Return the correlation function"""
  if npol is None: npol = ne
  mus = get_moments_ij(m_in/scale,n=npol,i=i,j=j,use_fortran=use_fortran)
  if np.sum(np.abs(mus.imag))>0.001:
#    print("WARNING, off diagonal has nonzero imaginary elements",np.sum(np.abs(mus.imag)))
    pass
  if x is None: xs = np.linspace(-1.0,1.0,ne,endpoint=True)*0.99 # energies
  else: xs = x/scale # use from input
  ysr,ysi = generate_profile_pair(mus.real,mus.imag,xs,kernel="jackson",
          use_fortran=use_fortran) # so it is the Green function
  ys = (ysr - 1j*ysi)/scale*np.pi
  return (scale*xs,ys)



def dm_vivj_energy(m_in,vi,vj,scale=10.,npol=None,ne=500,x=None):
  """Return the correlation function"""
  if npol is None: npol = ne
  mus = get_moments_vivj(m_in/scale,vi,vj,n=npol)
  if np.sum(np.abs(mus.imag))>0.001:
#    print("WARNING, off diagonal has nonzero imaginary elements",np.sum(np.abs(mus.imag)))
    pass
  if x is None:
    xs = np.linspace(-1.0,1.0,npol*10,endpoint=True)*0.95 # energies
    ysr,ysi = generate_profile_pair(mus.real,mus.imag,xs,kernel="jackson",
            use_fortran=use_fortran) # so it is the Green function
    ys = (ysr - 1j*ysi)/scale*np.pi
    xs = scale*xs # reescale
    ys = ys/scale # reescale
    return xs,ys
  else:
    # Evaluate the Chebyshev reconstruction directly at the caller's
    # requested frequencies instead of building a dense npol*10-point
    # grid (spanning the full [-0.95,0.95]*scale domain, regardless of
    # how few/localized the requested points are) and interpolating
    # down to `x` afterwards -- for the common case (x has a few
    # hundred/thousand points, npol in the thousands) this used to be
    # the dominant cost of the whole KPM dynamical correlator. Points
    # outside the KPM-valid domain [-0.95,0.95]*scale are set to 0,
    # matching the previous interp1d(...,fill_value=0.0) behavior.
    x = np.asarray(x,dtype=float)
    xr = x/scale
    valid = np.abs(xr)<=0.95
    yout = np.zeros(x.shape,dtype=np.complex128)
    if np.any(valid):
      ysr,ysi = generate_profile_pair(mus.real,mus.imag,xr[valid],
              kernel="jackson",use_fortran=use_fortran)
      yout[valid] = (ysr - 1j*ysi)/scale*np.pi/scale
    return x,yout








# ---------------------------------------------------------------------
# High-order delta-Chebyshev (HODC) reconstruction
#
# Yi, Massatt, Horning, Luskin, Pixley & Kaye, "A high-order regularized
# delta-Chebyshev method for computing spectral densities",
# arXiv:2512.03149 (2025).
#
# Ordinary KPM damps the moments, mu_k -> g_k mu_k, and reconstructs with
# the Chebyshev expansion of delta itself; the damping factors g_k are
# what make the result a *positive* kernel of width ~pi/p, and the price
# is that every such kernel is only second-order accurate: the
# reconstructed density is rho convolved with something whose error at a
# smooth point is O(p^-2).
#
# HODC keeps the moments untouched and changes what is reconstructed.
# delta(E-x) is replaced by the order-m rational regularization (Eq. 15)
#
#     K_eta(E,x) = -1/pi sum_l Im[ w_l / (E - x + eta*z_l) ] ,
#
# a sum of m complex-weighted Lorentzians of width eta, with poles
# z_l = x_l + i (Eq. 16). Writing F(zeta) = int rho(x)/(zeta-x) dx for the
# Cauchy transform of the density and expanding F(E+eta*z_l) about
# eta = 0, the l-sum multiplies the j-th derivative by sum_l w_l z_l^j, so
# imposing the Vandermonde moment-matching conditions
#
#     sum_l w_l z_l^j = delta_{j0} ,   j = 0 ... m-1                 (*)
#
# annihilates every term below j = m and leaves K_eta -> delta weakly at
# O(eta^m). m = 1 is w = 1, z = i, i.e. plain Lorentzian broadening,
# which is only O(eta) -- the Lorentzian's second moment diverges, so the
# usual "even kernel kills the linear term" argument does not apply and
# the tails dominate. (*) also kills the large-|u| tail of K_eta itself
# through the same cancellation, leaving eta^m/u^(m+1) rather than the
# Lorentzian's eta/u^2.
#
# Because delta itself never appears, the reconstruction is no longer
# "damped moments times T_k(E)": expanding K_eta(E,.) in Chebyshev
# polynomials (Eq. 17) gives energy-*dependent* coefficients nu_k(E,eta),
# and (Eq. 18)
#
#     rho^eta(E) = sum_{k<p} nu_k(E,eta) mu_k ,    mu_k = <b|T_k(H)|a>,
#
# with the *raw* moments -- the same ones ordinary KPM starts from. So
# this is pure post-processing: no extra Chebyshev recursion, no change
# to any DMRG/ED backend, and both kernels can be compared from a single
# moment run.
#
# The paper computes nu_k by evaluating (15) at Chebyshev nodes and
# applying a fast cosine transform. That is unnecessary here: the
# Chebyshev expansion of a single Cauchy kernel is classical,
#
#     1/(zeta-x) = (1/sqrt(zeta^2-1)) [ 1 + 2 sum_{k>=1} q^k T_k(x) ],
#     q = zeta - sqrt(zeta^2-1)   (the root with |q| < 1),
#
# so with zeta_l = E + eta*z_l the coefficients are closed-form,
#
#     nu_k(E,eta) = (2-delta_{k0}) * (-1/pi) * Im sum_l w_l q_l^k
#                                                 / sqrt(zeta_l^2-1) ,
#
# and are real, so complex moments (a != b) work exactly as they do for
# jackson. |q_l| = 1 - O(eta) makes q_l^k the whole convergence story:
# the truncation at p terms costs ~exp(-p*eta), which is why eta cannot
# be taken arbitrarily small at fixed p -- see hodc_default_eta.
#
# Caveat from the paper, reproduced here: for m > 2 the kernel is not
# positive, so the reconstructed density can go slightly negative near
# sharp features. That is the trade for the higher order, not a bug.
# ---------------------------------------------------------------------

HODC_MAX_ORDER = 8 # solving (*) in double precision degrades beyond this

# p*eta at which the default eta is placed. This is not arbitrary: it is
# the value that dmrgpy's own KPM conventions already imply. kpmdmrg's
# get_dynamical_correlator asks for n = (emax-emin)/delta * kpm_n_scale
# polynomials and rescales the spectrum by scale = 1/((emax-emin)*
# kpm_scale), so setting eta to the *requested* resolution delta gives
# n*eta*scale = kpm_n_scale/kpm_scale = 3/0.7 = 4.3 with the default
# kpm_n_scale/kpm_scale. Empirically that also sits in the flat minimum of
# the error-vs-eta curve *for accurate moments* (see examples/
# dynamical_correlator/hodc_VS_jackson_kernel): below ~2 the truncated
# Chebyshev series rings, above ~8 the O(eta^m) smoothing error takes
# over. The optimum moves up when the moments themselves are noisy --
# nu_k ~ q^k decays like exp(-k*eta), so eta is also the only damping
# HODC applies to the high-k moments, which are exactly the ones an MPS
# truncation gets wrong. That, not the value below, is the knob to reach
# for when an HODC spectrum from a modest-bond-dimension run looks worse
# than the Jackson one.
HODC_DEFAULT_P_ETA = 4.0


def hodc_default_eta(npol):
    """Default HODC regularization width for npol Chebyshev moments, in
    the rescaled energy units where the spectrum lives in [-1,1]."""
    return HODC_DEFAULT_P_ETA/float(npol)


def hodc_poles_weights(order=6):
    """Poles z_l = x_l + i and weights w_l of the order-m regularized
    delta (Eq. 15/16 of arXiv:2512.03149), solving the Vandermonde
    moment-matching system sum_l w_l z_l^j = delta_{j0}, j < m.

    The paper fixes Im z_l = 1 but leaves the real parts x_l unspecified.
    They are taken here as `order` equispaced points of unit spacing,
    symmetric about 0 -- the spacing matched to the (unit) pole height,
    which is the balanced choice: crowding the x_l together makes the
    Vandermonde system singular and the weights blow up, while spreading
    them out widens the kernel at fixed eta for no gain. The symmetry
    about 0 makes the pole set invariant under z -> -conj(z), hence
    w -> conj(w) by uniqueness, hence K_eta even in E-x."""
    order = int(order)
    if order<1 or order>HODC_MAX_ORDER:
        raise ValueError("hodc: order must be in 1..%d, got %s"%(
            HODC_MAX_ORDER,order))
    z = (np.arange(order) - (order-1)/2.0) + 1j
    V = np.vander(z,N=order,increasing=True).T # V[j,l] = z_l**j
    rhs = np.zeros(order,dtype=np.complex128) ; rhs[0] = 1.0
    return z,np.linalg.solve(V,rhs)


def _hodc_qa(xs,order,eta):
    """(q_l(E), a_l(E)) with a_l = w_l/sqrt(zeta_l^2-1), shape (m,len(xs))."""
    z,w = hodc_poles_weights(order)
    zeta = np.asarray(xs,dtype=np.complex128)[None,:] + eta*z[:,None]
    # sqrt(zeta-1)*sqrt(zeta+1), not sqrt(zeta**2-1): the factored form
    # avoids the cancellation in zeta**2-1 for |zeta| close to 1, which is
    # exactly where the KPM window's endpoints sit.
    s = np.sqrt(zeta-1.0)*np.sqrt(zeta+1.0)
    q = zeta - s
    q = np.where(np.abs(q)>1.0,zeta+s,q) # pick the root inside the disc
    return q,w[:,None]/s


def hodc_coefficients(xs,npol,order=6,eta=None):
    """nu_k(E,eta) of arXiv:2512.03149 Eq. (17), as a real (npol,len(xs))
    array, so that rho(E) = sum_k nu_k(E) mu_k.

    This materializes the full coefficient matrix and is meant for
    testing/inspection; generate_profile_hodc() evaluates the same sum
    without ever storing it."""
    xs = np.asarray(xs)
    npol = int(npol)
    if eta is None: eta = hodc_default_eta(npol)
    q,a = _hodc_qa(xs,order,float(eta))
    nus = np.zeros((npol,xs.size))
    qp = np.ones_like(q)
    for k in range(npol):
        nus[k] = -np.sum(a*qp,axis=0).imag/np.pi
        if k>0: nus[k] *= 2.0
        qp = qp*q
    return nus


def _hodc_horner(cs,q,a):
    """sum_k (2-delta_k0) cs[k] * (-1/pi) Im sum_l a_l q_l^k, for a real
    coefficient array cs. Since cs is real it commutes with Im, so the
    whole k-sum collapses to one Horner evaluation of
    P(q) = cs[0] + 2 sum_{k>=1} cs[k] q^k at the m poles -- O(m*npol)
    work per energy and no (npol,npts) temporary."""
    n = len(cs)
    b = np.zeros_like(q)
    for k in range(n-1,0,-1): b = b*q + 2.0*cs[k]
    b = b*q + cs[0]
    return -np.sum(a*b,axis=0).imag/np.pi


def generate_profile_hodc(mus,xs,order=6,eta=None):
    """Reconstruct a spectral density from raw Chebyshev moments with the
    high-order delta-Chebyshev kernel instead of a damped-moment KPM one.

    mus are the *undamped* moments <b|T_k(H)|a> (exactly what
    generate_profile receives), xs the rescaled energies in [-1,1], order
    the accuracy order m, and eta the regularization width in the same
    rescaled units (default hodc_default_eta(len(mus))). Returns a complex
    array, matching generate_profile's convention."""
    mus = np.asarray(mus)
    xs = np.asarray(xs)
    npol = len(mus)
    if eta is None: eta = hodc_default_eta(npol)
    eta = float(eta)
    if eta<=0.0: raise ValueError("hodc: eta must be positive, got %s"%eta)
    q,a = _hodc_qa(xs,order,eta)
    if np.iscomplexobj(mus): # nu_k is real, so real/imag parts decouple
        return (_hodc_horner(mus.real,q,a)
                + 1j*_hodc_horner(mus.imag,q,a))
    return _hodc_horner(np.asarray(mus,dtype=float),q,a) + 0.0j


def generate_profile(mus,xs,kernel="jackson",use_fortran=use_fortran,
        hodc_order=6,hodc_eta=None):
    """ Uses the Chebychev expansion to create a certain profile"""
    # HODC is not a damped-moment kernel -- it replaces the whole
    # reconstruction (energy-dependent coefficients, no 1/sqrt(1-x^2)
    # prefactor), so it short-circuits everything below, the fortran
    # path included. See generate_profile_hodc.
    if kernel=="hodc":
        return generate_profile_hodc(mus,xs,order=hodc_order,eta=hodc_eta)
    # initialize polynomials
  #  xs = np.array([0.])
    tm = np.zeros(xs.shape) +1.
    t = xs.copy()
    if kernel=="jackson": mus = jackson_kernel(mus)
    elif kernel=="lorentz": mus = lorentz_kernel(mus)
    elif kernel=="plain": pass # do nothing
    elif kernel is None: pass
    else: raise
    if use_fortran: # call the fortran routine
      ys = kpmf90.generate_profile(mus,xs) 
    else: # do a python loop
      ys = np.zeros(xs.shape,dtype=np.complex128) + mus[0] # first term
      # loop over all contributions
      for i in range(1,len(mus)):
        mu = mus[i]
        ys += 2.*mu*t # add contribution
        tp = 2.*xs*t - tm # chebychev recursion relation
        tm = t + 0.
        t = 0. + tp # next iteration
      ys = ys/np.sqrt(1.-xs*xs) # prefactor
    ys = ys/np.pi
    return ys



def generate_profile_pair(mus_r,mus_i,xs,kernel="jackson",use_fortran=use_fortran,
        hodc_order=6,hodc_eta=None):
    """Evaluate generate_profile for two real moment arrays (mus_r, mus_i)
    against the same xs grid in one pass. generate_profile's Chebyshev
    recursion (tp = 2*xs*t - tm) depends only on xs, not on the moments,
    so the callers that need mus.real and mus.imag reconstructed
    separately (dm_ij_energy/dm_vivj_energy -- kept apart rather than
    summed, since the final result combines them as ysr - 1j*ysi, not
    generate_profile(mus_r+1j*mus_i)) used to redo that recursion twice;
    this shares it across both instead."""
    if kernel=="hodc": # see generate_profile's own note
        return (generate_profile_hodc(mus_r,xs,order=hodc_order,eta=hodc_eta),
                generate_profile_hodc(mus_i,xs,order=hodc_order,eta=hodc_eta))
    if kernel=="jackson": mus_r = jackson_kernel(mus_r); mus_i = jackson_kernel(mus_i)
    elif kernel=="lorentz": mus_r = lorentz_kernel(mus_r); mus_i = lorentz_kernel(mus_i)
    elif kernel=="plain": pass # do nothing
    elif kernel is None: pass
    else: raise
    if use_fortran: # call the fortran routine (no fused variant available)
      ysr = kpmf90.generate_profile(mus_r,xs)
      ysi = kpmf90.generate_profile(mus_i,xs)
    else: # do a python loop, sharing the Chebyshev recursion
      ysr = np.zeros(xs.shape,dtype=np.complex128) + mus_r[0] # first term
      ysi = np.zeros(xs.shape,dtype=np.complex128) + mus_i[0] # first term
      tm = np.zeros(xs.shape) +1.
      t = xs.copy()
      for i in range(1,len(mus_r)):
        ysr += 2.*mus_r[i]*t # add contribution
        ysi += 2.*mus_i[i]*t # add contribution
        tp = 2.*xs*t - tm # chebychev recursion relation
        tm = t + 0.
        t = 0. + tp # next iteration
      denom = np.sqrt(1.-xs*xs) # prefactor
      ysr = ysr/denom
      ysi = ysi/denom
    return ysr/np.pi,ysi/np.pi


def generate_green_profile(mus,xs,kernel="jackson",use_fortran=use_fortran):
  """ Uses the Chebychev expansion to create a certain profile"""
  # initialize polynomials
#  xs = np.array([0.])
  tm = np.zeros(xs.shape) +1.
  t = xs.copy()
  ys = np.zeros(xs.shape,dtype=np.complex128) + mus[0]/2 # first term
  if kernel=="jackson": mus = jackson_kernel(mus)
  elif kernel=="lorentz": mus = lorentz_kernel(mus)
  else: raise
  if True:
    for i in range(1,len(mus)): # loop over mus
      ys += np.exp(1j*i*np.arccos(xs))*mus[i] # add contribution
    ys = ys/np.sqrt(1.-xs*xs)
    return 1j*2*ys/np.pi
#    if use_fortran: # call the fortran routine
#      ys = kpmf90.generate_profile(mus,xs) 
#    return ys








def dos(m_in,xs,ntries=20,n=200,scale=10.):
  """Return the density of states"""
  if scale is None: scale = 10.*np.max(np.abs(m_in.data)) # estimate of the value
  mus = random_trace(m_in/scale,ntries=ntries,n=n)
  ys = generate_profile(mus,xs/scale) # generate the DOS
  return ys # return the DOS 



def jackson_kernel(mus):
  """ Modify coeficient using the Jackson Kernel"""
  mo = mus.copy() # copy array
  n = len(mo)
  pn = np.pi/(n+1.) # factor
  for i in range(n):
    fac = ((n-i+1)*np.cos(pn*i)+np.sin(pn*i)/np.tan(pn))/(n+1)
    mo[i] *= fac
  return mo



def lorentz_kernel(mus):
  """ Modify coeficient using the Jackson Kernel"""
  mo = mus.copy() # copy array
  n = len(mo)
  pn = np.pi/(n+1.) # factor
  lamb = 3.
  for i in range(n):
    fac = np.sinh(lamb*(1.-i/n))/np.sinh(lamb)
    mo[i] *= fac
  return mo





def fejer_kernel(mus):
  """Default kernel"""
  n = len(mus)
  mo = mus.copy()
  for i in range(len(mus)):
    mo[i] *= (1.-float(i)/n) 
  return mo



def edge_dos(intra0,inter0,scale=4.,w=20,npol=300,ne=500,bulk=False,
                use_random=True,nrand=20):
  """Calculated the edge DOS using the KPM"""
  h = [[None for j in range(w)] for i in range(w)]
  intra = csc_matrix(intra0)
  inter = csc_matrix(inter0)
  for i in range(w): h[i][i] = intra
  for i in range(w-1): 
    h[i+1][i] = inter.H
    h[i][i+1] = inter
  h = bmat(h) # sparse hamiltonian
  ds = np.zeros(ne)
  dsb = np.zeros(ne)
  norb = intra0.shape[0] # orbitals ina cell
  for i in range(norb):
    (xs,ys) = ldos0d(h,i=i,scale=scale,npol=npol,ne=ne) 
    ds += ys # store
    if bulk:
      (xs,zs) = ldos0d(h,i=w*norb//2 + i,scale=scale,npol=npol,ne=ne) 
      dsb += zs # store
  if not bulk: return (xs,ds/w)
  else: return (xs,ds/w,dsb/w)


from .kpmextrapolate import extrapolate_moments,deconvolution



def reconstruct_chebyshev(mus,shift=0.,scale=1.0,
        x=np.linspace(-1.,1.,2000),kernel="jackson"):
    num_p = len(mus)
    xs2 = 0.99*np.linspace(-1.0,1.0,int(num_p*10),endpoint=False) # energies
    ys2 = generate_profile(mus,xs2,use_fortran=False,kernel=kernel) # generate the DOS
    xs2 += shift # add the shift
    xs2 *= scale # scale
    ys2 /= scale # scale
    if x is None: return xs2,ys2
    else:
      from scipy.interpolate import interp1d
      fr = interp1d(xs2, ys2.real,fill_value=0.0,bounds_error=False)
      fi = interp1d(xs2, ys2.imag,fill_value=0.0,bounds_error=False)
      return x,fr(x)+1j*fi(x)

