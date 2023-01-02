# methods to solve linear equations, purely written in Python
# they do not work properly yet

import numpy as np


# these functions are quite buggy, not functional at the current stage

def solve_Ab_cv(A, b,tol=1e-4,nmax=-1):
    raise # this is buggy
    x = b.copy()  # initialize x as the initial guess
    A = b.MBO.toMPO(A) # transform into an static MPO
    r = b - A*x  # compute the initial residual
    norm0 = b.norm() # input norm
    ii = 0
    print("here")
    while r.norm() > tol:  # continue until the residual is small enough
        alpha = r.dot(r) / r.dot(A*r)  # compute the correction factor
        x = x + alpha * r  # update the solution
        c = A*x
        norm1 = c.norm() # compute the norm of the solution
        x = x/norm1*norm0 # renormalize
        r = b - c  # update the residual
        ii +=1
        print(r.norm(),np.abs(alpha))
        if ii==nmax: break
    return x


def solve_Ab_cg(A, b,tol=1e-4,nmax=-1):
    raise # this is buggy
    def dot(x, y):
        return x.dot(y)

    norm0 = b.norm() # input norm
    A = b.MBO.toMPO(A) # transform into an static MPO
    x = b.copy()  # initialize x as the initial guess
    r = b - A*x  # compute the initial residual
    p = r.copy()  # initialize the search direction as the residual
    print("Starting")
    ii = 0
    while True:  # continue until the residual is small enough
        alpha = dot(r, r) / dot(p, A*p)  # compute the correction factor
        x = x + alpha * p  # update the solution
        r_new = r - alpha * (A*p)  # update the residual
        beta = dot(r_new, r_new) / dot(r, r)  # compute the search direction update factor
        p = r_new + beta * p  # update the search direction
        r = r_new  # update the residual
        norm1 = (A*x).norm() # input norm
        x = x*(norm0/norm1) # renormalize
        print(r.norm(),alpha)
        if r.norm() < tol:  # check if the residual is small enough
            break
        ii += 1
        if ii==nmax: break

    return x



def bicgstab(A, b, x0 = None, tol=1e-4, nmax=1e3):
    raise # this is buggy
    if x0 is None: x_0 = b.copy()
    else: x_0 = x0
    max_iter = nmax
    r = b - A * x_0
    r_tilde = r.copy()
    p = r.copy()
    p_tilde = r_tilde.copy()
    x = x_0.copy()
    
    r_dot_r_tilde = r.dot(r_tilde)
    r_norm = np.sqrt(r_dot_r_tilde)
    r_tol = tol * r_norm
    norm0 = b.norm() # initial norm
    num_iter = 0
    while r_norm > r_tol and num_iter < max_iter:
        num_iter += 1
        
        Ap = A * p
        alpha = r_dot_r_tilde / (p_tilde.dot(Ap))
        x = x + alpha * p
        r = r - alpha * Ap
        
        r_tilde = r_tilde - alpha * (A * p_tilde)
        r_dot_r_tilde_new = r.dot(r_tilde)
        
        beta = r_dot_r_tilde_new / r_dot_r_tilde
        p = r + beta * p
        p_tilde = r_tilde + beta * p_tilde
        
        r_dot_r_tilde = r_dot_r_tilde_new
        r_norm = np.sqrt(r_dot_r_tilde)
        norm1 = (A*x).norm() # output norm
        x = x*(norm0/norm1) # renormalize
    
    return x






solve_Ab = bicgstab

