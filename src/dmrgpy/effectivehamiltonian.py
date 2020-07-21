import numpy as np
import scipy.linalg as lg

#def get_effective_hamiltonian(self,**kwargs):
#    """Return an effective Hamiltonian"""
#    n = self.ns # number of sites 
#    pairs = []
#    for i in range(n): # loop
#        for j in range(n): # loop
#            pairs.append([i,j]) # store
#    cs = self.get_correlator(pairs=pairs,apply_hamiltonian=True,**kwargs)
#    norm = self.get_correlator(pairs=pairs,apply_hamiltonian=False,**kwargs)
#    h = np.zeros((n,n),dtype=np.complex) # create matrix with coefficients
#    b = np.zeros((n,n),dtype=np.complex) # create matrix with overlaps
#    for k in range(len(cs)): # loop
#        i,j = pairs[k] # get index
#        h[i,j] = cs[k] # store
#        b[i,j] = norm[k] # store
#    return h,b # return effective hamiltonian


def get_effective_hamiltonian_coefficients(self,mode="DMRG",n=4,
        tol=1e-3,operators=None):
    """Return an effective Hamiltonian"""
    if operators is None: return NotImplemented
    (es,ws) = self.get_excited_states(n=n,mode=mode) # return the excited states
    ops = dict() # different operators considered
    opm = dict() # representation of the different operators considered
    for i in range(len(operators)):  ops[i] = operators[i] # convert to dict
    for key in ops: opm[key] = get_representation(ws,ops[key])
    h = get_representation(ws,self.hamiltonian) # return the representation
    # now fit the matrix with the other ones
    opm[("Id")] = np.identity(len(es)) # identity
    coef = fit_matrix(h,opm) # fit the matrix with that dictionary
    del coef[("Id")] # remove identity
    return [coef[key] for key in coef] # return list



def get_effective_hamiltonian_couplings(self,mode="DMRG",n=4,method="single",
        tol=1e-3,operators=None):
    """Return an effective Hamiltonian"""
    (es,ws) = self.get_excited_states(n=n,mode=mode) # return the excited states
    ops = dict() # different operators considered
    op = get_projection_operators(self) # return the operators
    opm = dict() 
    if method=="single": # compute each operator
      for key in op: opm[key] = get_representation(ws,op[key])
    ops[("Id")] = np.identity(len(es)) # identity
    # use just by quadratic operators
    for key1 in op:   # loop
      for key2 in op:  # loop
          if method=="single":
              m = opm[key1]@opm[key2] # get the operator
          elif method=="full":
              o = op[key1]*op[key2] # get the operator
              m = get_representation(ws,o) # get the representation 
          else: raise
          if acceptable_matrix(m,ops): # if the matrix can be accepted
              ops[(key1,key2)] = m # store this matrix
    h = get_representation(ws,self.hamiltonian) # return the representation
    # now fit the matrix with the other ones
    coef = fit_matrix(h,ops) # fit the matrix with that dictionary
    del coef[("Id")] # remove identity
    return coef



def get_effective_hamiltonian(self,tol=1e-3,**kwargs):
    """Compute the effective Hamiltonian and return its latex form"""
    coef = get_effective_hamiltonian_couplings(self,**kwargs)
    return dict2latex(coef,tol=tol) # return the latex form



def dict2latex(d,tol=1e-3):
    """Transform the dictionary into a latex form"""
    cs = [d[key] for key in d] # coefficients
    cmax = np.round(np.max(np.abs(cs)),3) # maximum value
    out = "H = \n"+str(cmax)+" \[ \n" # output string
    for key in d: # loop
        c = np.round(d[key]/cmax,3) # round the number
        if np.abs(c)<tol: continue
        out += str(c) + "  " # normalize
        for k in key: out += k + "  "
        out += " + \n" # new line
    out += " \] \n" # last line
    return out






def get_projection_operators(self,mode="spin"):
    op = dict() # dictionary
    ns = self.ns//2 # number of sites
    for i in range(ns): 
        op["S^x_"+str(i+1)] = self.Sx[i]
        op["S^y_"+str(i+1)] = self.Sy[i]
        op["S^z_"+str(i+1)] = self.Sz[i]
    return op # return operators


def get_representation(ws,op):
    """Return the representation of a certain operator"""
    ne = len(ws)
    h = np.zeros((ne,ne),dtype=np.complex)
    for i in range(ne):
        for j in range(ne):
            h[i,j] = ws[i].overlap(op*ws[j])
    return h # return matrix


def acceptable_matrix(m,ops):
    """Check if it is ok to keep this matrix"""
    if np.sum(np.abs(m))<1e-7: return False
    v = matrix2vector(m)
    for key in ops: # loop over the other matrices
        o = ops[key] # get the matrix
        vo = matrix2vector(o) # convert to vector
        proj = v.dot(vo)/(np.sqrt(v.dot(v))*np.sqrt(vo.dot(vo)))
        if np.abs(proj)>0.98: 
      #      print("Skipping")
            return False
    return True



def matrix2vector(m):
    n = m.shape[0] # get the dimension
    v = np.zeros(2*n**2) # to a vector
    v[0:n**2] = m.reshape(n**2).real
    v[n**2:2*n**2] = m.reshape(n**2).imag
    return v



def fit_matrix(h,d,cutoff=1e-4):
    """Fit a matrix with a dictionary of matrices"""
    ms = np.array([d[key] for key in d]) # redefine as array
    n = len(ms) # number of matrices
    def f(v): # function to minimize
        diff = h.copy() # initialize
        rv = v[0:n] # real part
#        iv = v[n:2*n] # imaginary part
#        zv = rv+1j*iv # complex vector
        zv = rv
        for i in range(len(ms)): # loop over ms
            diff = diff - zv[i]*ms[i] # add this contribution
        error = np.mean(np.abs(diff)**2) # error
#        print("Error",error)
        return error
    from scipy.optimize import minimize
    x0 = np.random.random(n) # random guess
    sol = minimize(f,x0)
    x = sol.x # solution of the minimization
    x = x[0:n] #+ 1j*x[n:2*n] # redefine as complex
    out = dict()
    ii = 0
    h0 = 0.0
    for key in d: # loop over the operators 
        if np.abs(x[ii])>cutoff:
          out[key] = x[ii]
          h0 = h0 + x[ii]*d[key]
        ii += 1 # increase counter
#    print(np.round(h,2),"Original Hamiltonian")
#    print(np.round(h0,2),"Computed Hamiltonian")
    return out # return the coefficients








