"""The Kondo screening cloud of an Anderson impurity at the end of a
real-space tight-binding lead, as computed by Holzner, McCulloch,
Schollwoeck, von Delft and Heidrich-Meisner, PRB 80, 205114 (2009),
arXiv:0906.2933.

The model is the paper's eq. 2, a dot on site 0 and a lead on sites 1 to
L-1, so L counts every site, dot included:

    H = ed sum_s n_ds + B Sz_d + U n_dup n_ddn
        - t sum_{i=1}^{L-2} sum_s (c^dag_is c_i+1s + h.c.)
        - sqrt(2) t' sum_s (c^dag_1s d_s + h.c.),   Gamma = 2 t'^2/t

The paper writes the lead hopping sum up to L-1 but the total spin up to
site L-1; we follow the second, so the system has L sites and the sum rule
below closes at x = L-1.

The cloud is read off the integrated spin correlation, the paper's eq. 3,

    Sigma(x) = 1 + sum_{i=1}^x <S_d.S_i>/<S_d.S_d>

Sigma(0) = 1, and in a total-spin singlet <S_tot^2> = 0 gives
<S_d^2> = <(sum_i S_i)^2>, so sum_{i>=1} <S_d.S_i> = -<S_d^2> and
Sigma(L-1) = 0 exactly. That is why the paper uses an even L at half
filling: with L odd the ground state is a doublet and the sum rule fails.
In a field the paper uses the connected zz form,
(<Sz_d Sz_i> - <Sz_d><Sz_i>)/(<Sz_d Sz_d> - <Sz_d>^2), which has no sum
rule.

The screening length is the paper's eq. 5, the first distance at which a
fraction a of the moment is screened,

    xi_a(L) = min{x : Sigma(x) <= 1-a},

with a = 0.9 the paper's choice. It depends on L, since the sum rule forces
full screening at x = L-1, and it saturates once L is larger than the
cloud. The reference length is eq. 4, xi_K^0 = hbar v_F/T_K with the
Bethe-ansatz (Haldane) Kondo temperature of the Anderson model,

    T_K = sqrt(U Gamma) exp(-pi |ed| |ed+U|/(2 U Gamma)),

and v_F = 2t at half filling. The paper finds xi_a = p(a) xi_K^0/(hbar v_F)
= p(a)/T_K in units of t and of the lattice constant, with p(0.9) = 6.7,
for U/Gamma above roughly 2. The scaling collapses of the paper (Sigma
against x/xi_K, xi_a(L)/xi_a against L/xi_a) are fits done on these
numbers and are left to the user.

Every quantity is an expectation value, so mode="ED" cross-checks it on a
short chain. The paper does not use a Wilson chain or a modified hopping
to reach small T_K, it uses long uniform leads (L up to 500, SU(2)
DMRG), so neither is offered here; the cloud then fits only for U/Gamma
of a few.
"""

import warnings

import numpy as np


def anderson_lead(L, ed=-0.5, U=1.0, Gamma=0.2, t=1.0, B=0.0,
                  chain="native", half_filling=False, **kwargs):
    """Anderson impurity on site 0 of a lead of L-1 sites (L sites in
    total), the paper's eq. 2, as a spinful fermion chain.

    Gamma = 2 t'^2/t is the hybridization, the dot couples to site 1 by
    sqrt(2) t' = sqrt(Gamma t). B is a field B Sz on the dot only.

    chain is "native" (Spinful_Fermionic_Chain_Native, one site of
    dimension 4 per orbital) or "interleaved" (Spinful_Fermionic_Chain,
    two spinless sites per orbital). Native is the default: at L=40 and
    equal maxm it gave a lower energy and a sum-rule error 16 times
    smaller on "python", at about twice the time.

    The paper works at half filling, N = L electrons, for every ed. At
    ed = -U/2 particle-hole symmetry puts the ground state there by
    itself. Away from it, pass half_filling=True, which sets the conserved
    sector Nf = L. Take L even: at half filling the ground state is then
    a singlet. The layout is stored as fc.impurity and fc.channels, which
    get_kondo_cloud and get_kondo_screening read. Extra keywords go to the
    chain."""
    from . import fermionchain
    if L < 2:
        raise ValueError("the lead needs at least one site, L >= 2")
    if Gamma < 0.0:
        raise ValueError("Gamma must be non-negative")
    if L % 2 == 1:
        warnings.warn("L = %d sites is odd: at half filling the ground state "
            "then has an odd number of electrons and is a doublet, not a "
            "singlet, and Sigma(L-1) = 0 does not hold; take L even" % L)
    cls = {"interleaved": fermionchain.Spinful_Fermionic_Chain,
           "native": fermionchain.Spinful_Fermionic_Chain_Native}
    if chain not in cls:
        raise ValueError("chain must be 'interleaved' or 'native'")
    fc = cls[chain](L, **kwargs)
    hop = lambda i, j: (fc.Cdagup[i]*fc.Cup[j] + fc.Cdagdn[i]*fc.Cdn[j]
                        + fc.Cdagup[j]*fc.Cup[i] + fc.Cdagdn[j]*fc.Cdn[i])
    h = ed*(fc.Nup[0] + fc.Ndn[0]) + U*fc.Nup[0]*fc.Ndn[0]
    if B != 0.0:
        h = h + B*fc.Sz[0]
    h = h - np.sqrt(Gamma*t)*hop(0, 1)
    for i in range(1, L-1):
        h = h - t*hop(i, i+1)
    fc.set_hamiltonian(h)
    if half_filling:
        fc.set_conserved_sector(Nf=L)
    fc.impurity = [0]
    fc.channels = [list(range(1, L))]
    return fc


def kondo_temperature(U, Gamma, ed=None):
    """T_K = sqrt(U Gamma) exp(-pi |ed||ed+U|/(2 U Gamma)), the Kondo
    temperature that the paper's eq. 4 uses (ed = -U/2 if not given)"""
    if ed is None:
        ed = -U/2.0
    return np.sqrt(U*Gamma)*np.exp(-np.pi*abs(ed)*abs(ed+U)/(2.0*U*Gamma))


def kondo_length(U, Gamma, ed=None, t=1.0):
    """xi_K^0 = hbar v_F/T_K with v_F = 2t, the paper's eq. 4, in units of
    the lattice constant"""
    return 2.0*t/kondo_temperature(U, Gamma, ed=ed)


def screening_length(sigma, a=0.9):
    """xi_a = min{x : Sigma(x) <= 1-a}, the paper's eq. 5, with sigma[x]
    = Sigma(x). Returns nan when Sigma never drops to 1-a."""
    below = np.nonzero(np.asarray(sigma) <= 1.0 - a)[0]
    return int(below[0]) if len(below) > 0 else np.nan


def _dot(self, i, j):
    """S_i.S_j with real coefficients, Sz Sz + (S+ S- + S- S+)/2"""
    sp = lambda k: self.Cdagup[k]*self.Cdn[k]
    sm = lambda k: self.Cdagdn[k]*self.Cup[k]
    return self.Sz[i]*self.Sz[j] + 0.5*(sp(i)*sm(j) + sm(i)*sp(j))


def kondo_cloud(self, imp=None, lead=None, a=0.9, connected=False,
                mode="DMRG", sum_rule_tol=1e-3):
    """Spin correlations between the impurity and the lead, the integrated
    correlation Sigma(x) and the screening length xi_a, see the module
    docstring.

    imp is the impurity orbital and lead the lead orbitals ordered by
    distance; both default to the layout anderson_lead stores. a is a
    number or a list of numbers. connected=True takes the zz connected form
    the paper uses in a field. Without it, a warning says when
    |Sigma(L-1)| > sum_rule_tol, since then the state is not a singlet.
    Returns a dict with
      x      distances 0..len(lead)
      SdSi   <S_d.S_i> over the lead (connected zz form if connected)
      SdSd   <S_d.S_d> (or <Sz_d^2> - <Sz_d>^2)
      Sigma  Sigma(x), Sigma[0] = 1
      xi     xi_a (one number, or a dict a -> xi_a)
      sum_rule  Sigma(len(lead)), zero in a singlet when not connected
      mu     <Sz_d^2> - <Sz_d>^2, the moment of the dot
      nd     <n_d>"""
    if imp is None:
        imp = getattr(self, "impurity", [0])[0]
    if lead is None:
        ch = getattr(self, "channels", None)
        lead = ch[0] if ch is not None else [i for i in range(len(self.Sz))
                                             if i != imp]
    lead = list(lead)
    vev = lambda op: np.real(self.vev(op, mode=mode))
    szd = vev(self.Sz[imp])
    sz2 = vev(self.Sz[imp]*self.Sz[imp])
    mu = sz2 - szd**2
    if connected:
        norm = mu
        corr = np.array([vev(self.Sz[imp]*self.Sz[i]) - szd*vev(self.Sz[i])
                         for i in lead])
    else:
        norm = vev(_dot(self, imp, imp))
        corr = np.array([vev(_dot(self, imp, i)) for i in lead])
    if norm <= 0.0:
        raise ValueError("the impurity carries no moment, <S_d^2> = 0, so "
                         "Sigma(x) is not defined")
    sigma = 1.0 + np.concatenate([[0.0], np.cumsum(corr)])/norm
    if np.ndim(a) == 0:
        xi = screening_length(sigma, a)
    else:
        xi = {ai: screening_length(sigma, ai) for ai in a}
    if not connected and abs(sigma[-1]) > sum_rule_tol:
        warnings.warn("Sigma(L-1) = %.2e, not 0: the state is not a total-"
            "spin singlet (odd L, a filling away from N = L, or a maxm too "
            "small), so the sum rule does not close and xi_a is not the "
            "paper's" % sigma[-1])
    nd = vev(self.Nup[imp] + self.Ndn[imp])
    return {"x": np.arange(len(lead)+1), "SdSi": corr, "SdSd": norm,
            "Sigma": sigma, "xi": xi, "sum_rule": sigma[-1], "mu": mu,
            "nd": nd}
