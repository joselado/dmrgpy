from . import operatornames
import numpy as np


# On itensor_version="python" the profiles and the two-point loops below
# go through one sweep of the state (pyitensor/chain.py's sweep_measure,
# the algorithm of ITensorMPS.jl's expect and correlation_matrix) instead
# of one MPO and one contraction of the whole chain per vev. Set to False
# to take the per-term vev on every backend: that is the reference the
# sweep is checked against (tests/test_pyitensor_expect_sweep.py).
ONE_SWEEP = True


def _route(self, kwargs):
    """(session, state) for the one-sweep route, or None when the per-term
    vev answers: another backend, mode="ED" (asked or forced), a chain too
    short for DMRG, or a keyword other than mode= and wf= (npow=, T=, or
    one vev() refuses)."""
    if not ONE_SWEEP or getattr(self, "itensor_version", None) != "python":
        return None
    if set(kwargs) - {"mode", "wf"}:
        return None
    if getattr(self, "ns", 0) < 2:
        return None  # get_mode() falls back to ED, and prints it
    if self.get_mode(mode=kwargs.get("mode", "DMRG")) != "DMRG":
        return None
    session = getattr(self, "_session", None)
    if not hasattr(session, "sweep_measure"):
        return None
    wf = kwargs.get("wf")
    if wf is None:
        wf = self.get_gs()  # as multi_vev does
    return session, wf.cpp_handle


def _pieces(op):
    """op as a sum of products on one or two sites: [(c, groups)], groups
    a tuple of (0-based site, names) in the written order, names the tuple
    of one-site names multiplied there. An empty groups is a constant.
    None when a term touches three sites or more, or comes back to a site
    it has left. The names are those the operator was built with, before
    the Jordan-Wigner transform of to_terms(), since the sweep threads the
    string as AutoMPO does."""
    from . import multioperator
    if type(op) != multioperator.MultiOperator:
        return None
    out = []
    for c, factors in op.to_terms(jordan_wigner_transform=False):
        groups = []
        for name, site in factors:
            if name == "Id":
                continue
            s = site - 1
            if groups and groups[-1][0] == s:
                groups[-1][1].append(name)
            elif len(groups) == 2 or any(g[0] == s for g in groups):
                return None
            else:
                groups.append((s, [name]))
        out.append((c, tuple((s, tuple(names)) for s, names in groups)))
    return out


def sweep_vevs(self, ops, **kwargs):
    """<op> for every operator of `ops`, normalized as vev() normalizes,
    from one sweep of the state on itensor_version="python". None when
    the per-term vev has to answer: see _route, plus an operator with a
    term on three sites or more, a fermion-odd one-site term or pair (its
    string runs to the end of the chain), and, in a conserved sector, a
    term that does not conserve the sector on its own (vev() then raises,
    or answers a sum such as Sx Sx + Sy Sy that conserves it as a sum)."""
    ops = list(ops)
    if not ops:
        return None
    route = _route(self, kwargs)
    if route is None:
        return None
    session, wf = route
    split = [_pieces(op) for op in ops]
    if any(p is None for p in split):
        return None
    from .pyitensor.sites.base import is_fermionic
    def odd(names):
        return sum(1 for x in names if is_fermionic(x)) % 2 == 1
    ones, twos = {}, {}
    for pieces in split:
        for c, groups in pieces:
            factors = [(x, s+1) for s, names in groups for x in names]
            if not session.conserves_sector(factors):
                return None
            if len(groups) == 1:
                (s, p), = groups
                if odd(p):
                    return None
                ones.setdefault(p, {})[s] = None
            elif len(groups) == 2:
                (i, p), (j, q) = groups
                if odd(p) != odd(q):
                    return None
                twos.setdefault((p, q), {})[(i, j)] = None
    ones = [(p, list(sites)) for p, sites in ones.items()]
    twos = [(p, q, list(pairs)) for (p, q), pairs in twos.items()]
    v1, v2 = session.sweep_measure(wf, ones=ones, twos=twos)
    one = {(p, s): v for (p, sites), vals in zip(ones, v1)
           for s, v in zip(sites, vals)}
    two = {(p, q, i, j): v for (p, q, pairs), vals in zip(twos, v2)
           for (i, j), v in zip(pairs, vals)}
    norm = session.overlap(wf, wf).real
    out = []
    for pieces in split:
        total = 0j
        for c, groups in pieces:
            if len(groups) == 0:
                total += c
            elif len(groups) == 1:
                (s, p), = groups
                total += c*one[(p, s)]/norm
            else:
                (i, p), (j, q) = groups
                total += c*two[(p, q, i, j)]/norm
        out.append(total)
    return np.array(out)


def _ladder_ss(self, pairs):
    """S_i.S_j as Sz Sz + (S+ S- + S- S+)/2 for every pair, the form whose
    every term conserves Sz, so that the sweep also serves an Sz sector.
    None unless Sx, Sy and Sz are the sites' own one-site operators and
    every site of the session defines S+ and S-."""
    sites = getattr(getattr(self, "_session", None), "sites", None)
    if sites is None:
        return None
    for k in range(self.ns):
        for name, ops in (("Sx", self.Sx), ("Sy", self.Sy), ("Sz", self.Sz)):
            if _pieces(ops[k]) != [(1.0, ((k, (name,)),))]:
                return None
        try:
            sites.site_type(k+1).matrix("S+")
            sites.site_type(k+1).matrix("S-")
        except ValueError:
            return None
    op = self.get_operator
    return [op("Sz", i)*op("Sz", j)
            + 0.5*(op("S+", i)*op("S-", j) + op("S-", i)*op("S+", j))
            for (i, j) in pairs]


def get_correlator(self,pairs=[[]],name="SS",**kwargs):
    """Compute a certain static correlator in the spin chain"""
    if name=="SS":
        def getop(i,j):
            m = self.Sx[i]*self.Sx[j]
            m = m + self.Sy[i]*self.Sy[j]
            m = m + self.Sz[i]*self.Sz[j]
            return m
    else:
      namei,namej = operatornames.recognize(name) # return that one
      opi = operatornames.name2MO(namei,self)
      opj = operatornames.name2MO(namej,self)
      def getop(i,j):
          return opi[i]*opj[j]
    ops = [getop(i,j) for (i,j) in pairs]
    ladder = None
    if name=="SS" and ONE_SWEEP and getattr(self,"itensor_version",None)=="python":
        ladder = _ladder_ss(self,pairs)
    out = sweep_vevs(self,ops if ladder is None else ladder,**kwargs)
    if out is not None: return out
    return np.array([self.vev(o,**kwargs) for o in ops])
