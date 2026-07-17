import numpy as np
from . import operatornames
from . import multioperator

def dynamical_correlator(self,es=np.linspace(0.,10.0,100),
        delta=1e-1,name="XX",i=0,j=0):
    """
    Compute the dynamical correlator using CVM method in DMRG
    """
    if not self.computed_gs: self.get_gs() # compute ground state
    out = [] # empty list
    for e in es: # loop over energies
        print("CVM in E = ",e)
        o = cvm_dmrg(self,name=name,i=i,j=j,delta=delta,e=e)
        out.append(o) # store
    out = np.array(out)
#    from .inference import points2function
#    (es,out) = points2function(es,out)
    return (es,out) # return result










def cvm_dmrg(self,name="XX",delta=1e-1,e=0.0,**kwargs):
    """
    Return the dynamical correlator for a single energy/frequency using
    the Correction Vector Method, computed entirely in Python (see
    cvm_correction_vector below) rather than via the in-process pybind11
    extension's Chain::cvm_dynamical_correlator/bicstab (which is a hand-
    rolled, unguarded BiCGSTAB on the non-Hermitian z-H operator, prone to
    breakdown/blow-up for small eta -- see cvm_correction_vector's
    docstring). Parameter names differ from the underlying convention
    (omega, eta, energy) to match this file's existing naming: e -> omega,
    delta -> eta, self.e0 -> energy.
    """
    name = operatornames.str2MO(self,name,**kwargs)
    A = name[0]
    B = name[1]
    return cvm_correction_vector(self,A,B,e,delta,
            tol=self.cvm_tol,max_it=int(self.cvm_nit))


def cvm_correction_vector(self,A,B,omega,eta,tol=1e-5,max_it=1000):
    """
    Correction Vector Method (Ramasesha; Kuhner & White 1999):
    -Im<GS|A (omega+E0+i*eta-H)^{-1} B|GS>/pi, computed by solving the
    Hermitian, positive-definite system

        [(H-omega-E0)^2 + eta^2] xc = -eta * B|GS>

    via conjugate gradient, then recovering the correction vector as

        x = i*xc + (H-omega-E0)/eta * xc = (omega+E0+i*eta-H)^{-1} B|GS>

    This is the same algorithm the ED backend already uses
    (edtk/timedependent.py::solve_cv), just carried out with MPS algebra
    (StaticOperator/MPS arithmetic) instead of dense matrices. Because
    (H-omega-E0)^2+eta^2 is Hermitian positive-definite whenever H is
    Hermitian, *exact* CG would converge with a monotonically non-increasing
    residual and no breakdown mode -- unlike the previous direct solve of
    the non-Hermitian (z-H) system via a hand-rolled BiCGSTAB
    (Chain::cvm_dynamical_correlator/bicstab in chain_session.h), which had
    no protection against BiCGSTAB's well-known near-singular breakdown and
    could blow up. In practice every MPS here is truncated
    (maxdim=cvm_maxm) at every step, so the residual can still wander
    non-monotonically (confirmed directly: traced iteration-by-iteration on
    a 14-site chain) -- CG's breakdown-free guarantee is about avoiding
    BiCGSTAB's *divide-by-near-zero* failure mode, not about guaranteeing
    monotonic convergence once every intermediate vector is lossily
    compressed. To stay robust against that, the correction vector
    achieving the lowest residual seen so far is tracked and returned
    (rather than trusting whatever the loop ends on at max_it or an
    increased residual). Implemented purely with already-exposed Python
    primitives (self.toMPO()/StaticOperator for a build-once/apply-many
    MPO, MPS +/-/scalar-* and .dot() for the rest), so no new C++/pybind11
    bindings and no recompilation are needed to tune this further -- unlike
    the KPM/TDVP paths, which run their inner loop in C++.
    """
    wf0 = self.get_gs() # ground state (cheap/cached; also sets self.e0)
    self._session.set_sweep_params(self.cvm_maxm,self.nsweeps,self.cutoff,self.noise)
    self._session.set_verbose(self.verbose)
    self._session.set_mpomaxm(max(self.cvm_maxm,self.mpomaxm))
    Hshift = self.toMPO(self.hamiltonian-(self.e0+omega)) # (H-omega-E0), built once
    def applyA(v): return Hshift*(Hshift*v) + (eta*eta)*v # (H-omega-E0)^2+eta^2
    b = (-eta)*(B*wf0) # -eta*B|GS>
    xc = b # initial guess, matches the old bicstab's own x=b start
    r = b - applyA(xc)
    p = r
    rs_old = r.dot(r).real # ||r||^2, real since A is Hermitian PD
    best_xc,best_res = xc,np.sqrt(abs(rs_old))
    for k in range(max_it):
        Ap = applyA(p)
        alpha = rs_old/p.dot(Ap).real # real: <p|A|p> is real for Hermitian A
        xc = xc + alpha*p
        r = r - alpha*Ap
        rs_new = r.dot(r).real
        res = np.sqrt(abs(rs_new))
        if res<best_res: best_xc,best_res = xc,res # guard against truncation-driven non-monotonicity
        if res<=tol: break
        p = r + (rs_new/rs_old)*p
        rs_old = rs_new
    x = 1j*best_xc + (Hshift*best_xc)*(1./eta) # full correction vector
    G = wf0.dot(A*x) # <GS|A|x>
    return -G.imag/np.pi




def dynamical_correlator_analytic_continuation(self,name=None,
        delta=1e-1,es=np.linspace(0.,5.0,300)):
    """
    Compute the dynamical correlator using analytic continuation
    """
    A,B = name[0],name[1]
    wf = self.get_gs() # get the ground state
    wfa = A.get_dagger()*wf # apply A to the GS
    wfb = B*wf # apply B to the GS
    e0 = self.gs_energy() # ground state energy
    Hp = self.hamiltonian - e0
    def f(e): # function to compute
        wfi = self.applyinverse(-self.hamiltonian+(e0+e),wfa)
        return wfb.dot(wfi) # return result
#    return es,-np.array([f(e+1j*delta*10) for e in es]).imag*2/np.pi # brute force
    from .analyticcontinuation import imag2real
    xz = es*1j
    xz = np.linspace(delta*10,10.,100)*1j
    xz = np.concatenate([-xz,xz])
    xz = np.linspace(min(es),max(es),20) + delta*40*1j
#    xz = [np.random.random()-.5+1j*np.random.random()+0.5j for i in range(40)]
#    xz = 40.*np.array(xz)
    outz = np.array([f(z) for z in xz]) # complex axis
    esz,out = imag2real(xz,outz,x=es+1j*delta)
    out = -out.imag*2/np.pi
    return es,out



from .nonhermitian.dynamics import dynamical_correlator_cvm_explicit

#def dynamical_correlator_cvm_explicit(self,name=None,
#        delta=1e-1,es=np.linspace(0.,5.0,300)):
#    """
#    Compute the dynamical correlator using analytic continuation
#    """
#    ### So far this just works for onsite correlators
#    A,B = name[0],name[1]
#    if not self.is_zero_operator(A.get_dagger()-B): 
#        print("Only implemented for A^\dagger=B")
#        raise
#    wf = self.get_gs() # get the ground state
#    wfa = A.get_dagger()*wf # apply A to the GS
#    wfb = B*wf # apply B to the GS
#    e0 = self.gs_energy() # ground state energy
#    Hp = self.hamiltonian - e0
#    def f(e,delta): # function to compute
#        wfi = self.applyinverse(-self.hamiltonian+(e0+e+1j*delta),wfa)
#        return wfb.dot(wfi) # return result
#    from .analyticcontinuation import imag2real
#    outz = np.array([f(z,delta) - f(z,-delta) for z in es]) # complex axis
#    return es,1j*outz/np.pi





