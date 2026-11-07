"""The ED backend of mixedchain.Mixed_Spin_Fermion_Chain
(edtk/mixededchain.py): spin sites and spinful fermion locations in one
tensor-product Hilbert space, C/Cdag with a Jordan-Wigner string over the
fermion modes only.

What is pinned here:
- all-spin sitesin gives exactly Spin_Chain's ED spectrum, all-"F"
  sitesin exactly Spinful_Fermionic_Chain's ED spectrum;
- two genuinely mixed models, a spin 1/2 between two fermion locations
  (Kondo-like, hopping across the spin) and a spin 1 between fermion
  locations (hopping across the spin, Hubbard U, single-ion anisotropy),
  agree with the dense oracle of _dense_jw_reference.py to roundoff;
- the same two models agree between ED, "python" and v3 DMRG: energy,
  one-site vev, a spin-spin and a hopping correlator across the spin,
  a KPM dynamical correlator and <Sx(t)> of the spin after a quench;
- the ns<3 fallback of v3 and the ns<2 fallback of "python" answer with
  the ED number;
- conserved sectors: ED targets total Sz and Nf on a chain pinned to
  mode="ED", a sector the DMRG sites cannot carry is refused otherwise
  (submode="SECTOR" included), and the sectors the DMRG sites can carry
  agree with ED.
"""
import numpy as np
import pytest

from dmrgpy import mixedchain, spinchain, fermionchain

from _dense_jw_reference import (SX, SY, SZ, spin_one, embed, fermion_ops,
        gs as _gs)

BACKENDS = ("python", 3)
ED_TOL = 1e-10
E_TOL = 1e-6   # DMRG energy against ED
V_TOL = 1e-5   # DMRG vev and correlators against ED
QUENCH_TOL = 1e-6 # DMRG real-time evolution against ED


def _mc(labels, v="python", **kw):
    mc = mixedchain.Mixed_Spin_Fermion_Chain(labels, itensor_version=v, **kw)
    mc.maxm = 60
    mc.nsweeps = 20
    return mc


# ---------------------------------------------------------------------
# Degenerate sitesin: identical spectra to the dedicated ED backends
# ---------------------------------------------------------------------

def test_all_spin_ed_equals_spin_chain_ed():
    labels = ["1/2", "1", "1/2", "1"]
    def ham(c):
        h = 0.7*c.SS(0, 1) + 1.3*c.SS(1, 2) + 0.5*c.SS(2, 3)
        h = h + 0.4*c.Sz[0] - 0.9*c.Sz[1] + 0.2*c.Sx[3]
        return h + 0.3*c.Sz[1]*c.Sz[1] - 0.2*c.Sz[3]*c.Sz[3]
    sc = spinchain.Spin_Chain(labels)
    sc.set_hamiltonian(ham(sc))
    mc = _mc(labels)
    mc.set_hamiltonian(ham(mc))
    es_sc = np.array(sc.get_excited(n=8, mode="ED"))
    es_mc = np.array(mc.get_excited(n=8, mode="ED"))
    assert np.max(np.abs(es_sc - es_mc)) < ED_TOL
    for op in ("Sz", "Sx"):
        a = sc.vev(getattr(sc, op)[1]*getattr(sc, op)[2], mode="ED").real
        b = mc.vev(getattr(mc, op)[1]*getattr(mc, op)[2], mode="ED").real
        assert b == pytest.approx(a, abs=ED_TOL)
    kw = dict(es=np.linspace(-1., 3., 60), delta=0.1, submode="KPM", mode="ED")
    _, ya = sc.get_dynamical_correlator(name=(sc.Sz[1], sc.Sz[1]), **kw)
    _, yb = mc.get_dynamical_correlator(name=(mc.Sz[1], mc.Sz[1]), **kw)
    assert np.max(np.abs(ya - yb)) < 1e-9


def test_all_fermion_ed_equals_spinful_fermionic_chain_ed():
    n = 3
    def ham(c):
        h = 0
        for i in range(n-1):
            t = 0.6 + 0.3*i
            h = h + t*(c.Cdagup[i]*c.Cup[i+1] + c.Cdagdn[i]*c.Cdn[i+1])
        for i in range(n):
            h = h + 0.3*c.Delta[i] # onsite pairing, breaks Nf
        h = h + h.get_dagger()
        for i in range(n):
            h = h + (0.2*i - 0.1)*c.Ntot[i] + 1.1*c.Nup[i]*c.Ndn[i]
        return h + 0.25*c.Sx[0] # spin flip, breaks Sz
    fc = fermionchain.Spinful_Fermionic_Chain(n)
    fc.set_hamiltonian(ham(fc))
    mc = _mc(["F"]*n)
    mc.set_hamiltonian(ham(mc))
    es_fc = np.array(fc.get_excited(n=8, mode="ED"))
    es_mc = np.array(mc.get_excited(n=8, mode="ED"))
    assert np.max(np.abs(es_fc - es_mc)) < ED_TOL
    for op in (lambda c: c.Cdagup[0]*c.Cup[2], lambda c: c.Delta[1],
               lambda c: c.Sz[0]*c.Sz[2]):
        a = fc.vev(op(fc), mode="ED")
        b = mc.vev(op(mc), mode="ED")
        assert abs(a - b) < ED_TOL


# ---------------------------------------------------------------------
# Two genuinely mixed models, and their dense oracles
# ---------------------------------------------------------------------

KONDO = dict(K1=0.5, K2=0.35, t=0.4, muA=0.15, muB=-0.2, hf=0.13)
SPIN1 = dict(J1=0.8, J2=-0.4, D=0.6, t02=0.5, t23=0.7, U=1.5, mu=-0.6,
             hf=0.11)


def _kondo_mc(mc):
    """F - 1/2 - F, exchange to both fermion locations, hopping across"""
    p = KONDO
    h = (p["K1"]*mc.SS(0, 1) + p["K2"]*mc.SS(1, 2)
         + p["muA"]*mc.Ntot[0] + p["muB"]*mc.Ntot[2])
    hop = p["t"]*(mc.Cdagup[0]*mc.Cup[2] + mc.Cdagdn[0]*mc.Cdn[2])
    h = h + hop + hop.get_dagger()
    return h + p["hf"]*(mc.Sz[0] + mc.Sz[1] + mc.Sz[2])


def _spin1_mc(mc):
    """F - 1 - F - F, the spin 1 has a single-ion anisotropy D and sits
    between the first two fermion locations, which hop across it"""
    p = SPIN1
    h = p["J1"]*mc.SS(0, 1) + p["J2"]*mc.SS(1, 2) + p["D"]*mc.Sz[1]*mc.Sz[1]
    hop = p["t02"]*(mc.Cdagup[0]*mc.Cup[2] + mc.Cdagdn[0]*mc.Cdn[2])
    hop = hop + p["t23"]*(mc.Cdagup[2]*mc.Cup[3] + mc.Cdagdn[2]*mc.Cdn[3])
    h = h + hop + hop.get_dagger()
    for i in (0, 2, 3):
        h = h + p["U"]*mc.Nup[i]*mc.Ndn[i] + p["mu"]*mc.Ntot[i]
    return h + p["hf"]*(mc.Sz[0] + mc.Sz[1] + mc.Sz[2] + mc.Sz[3])


MODELS = {"kondo": (["F", "1/2", "F"], _kondo_mc),
          "spin1": (["F", "1", "F", "F"], _spin1_mc)}


def _oracle(model):
    """Dense Hamiltonian and observables of a model, built from
    _dense_jw_reference only. Locations: a fermion location is two modes
    (up, dn), a spin location one site."""
    labels, _ = MODELS[model]
    mask, dims, where = [], [], []
    for l in labels:
        if l == "F":
            where.append((len(mask), len(mask)+1))
            mask += [True, True]; dims += [2, 2]
        else:
            where.append(len(mask))
            mask.append(False); dims.append(2 if l == "1/2" else 3)
    ns = len(mask)
    C, Cd, N = fermion_ops(mask, dims=dims)
    S = {}
    for loc, l in enumerate(labels):
        if l == "F":
            u, d = where[loc]
            S[loc] = (0.5*(Cd[u]@C[d]) + 0.5*(Cd[d]@C[u]),
                      -0.5j*(Cd[u]@C[d]) + 0.5j*(Cd[d]@C[u]),
                      0.5*N[u] - 0.5*N[d])
        else:
            mats = (SX, SY, SZ) if l == "1/2" else spin_one()
            S[loc] = tuple(embed(m, where[loc], ns, dims=dims) for m in mats)
    SS = lambda i, j: sum(S[i][k]@S[j][k] for k in range(3))
    nt = lambda loc: N[where[loc][0]] + N[where[loc][1]]
    def hop(t, a, b):
        m = t*(Cd[where[a][0]]@C[where[b][0]] + Cd[where[a][1]]@C[where[b][1]])
        return m + m.conj().T
    if model == "kondo":
        p = KONDO
        H = (p["K1"]*SS(0, 1) + p["K2"]*SS(1, 2) + p["muA"]*nt(0)
             + p["muB"]*nt(2) + hop(p["t"], 0, 2)
             + p["hf"]*(S[0][2] + S[1][2] + S[2][2]))
    else:
        p = SPIN1
        H = (p["J1"]*SS(0, 1) + p["J2"]*SS(1, 2) + p["D"]*S[1][2]@S[1][2]
             + hop(p["t02"], 0, 2) + hop(p["t23"], 2, 3)
             + p["hf"]*sum(S[loc][2] for loc in range(4)))
        for loc in (0, 2, 3):
            u, d = where[loc]
            H = H + p["U"]*N[u]@N[d] + p["mu"]*nt(loc)
    nf = sum(nt(loc) for loc in range(len(labels)) if labels[loc] == "F")
    sz2 = 2*sum(S[loc][2] for loc in range(len(labels)))
    obs = {"Sz1": S[1][2], "N0": nt(0), "SzSz01": S[0][2]@S[1][2],
           "hop02": Cd[where[0][0]]@C[where[2][0]],
           "Nf": nf, "Sz2": sz2}
    return H, obs


def _observables(mc, **kw):
    return {"Sz1": mc.vev(mc.Sz[1], **kw), "N0": mc.vev(mc.Ntot[0], **kw),
            "SzSz01": mc.vev(mc.Sz[0]*mc.Sz[1], **kw),
            "hop02": mc.vev(mc.Cdagup[0]*mc.Cup[2], **kw)}


@pytest.mark.parametrize("model", sorted(MODELS))
def test_mixed_ed_matches_dense_oracle(model):
    labels, ham = MODELS[model]
    mc = _mc(labels)
    mc.set_hamiltonian(ham(mc))
    H, obs = _oracle(model)
    es_ref = np.linalg.eigvalsh(H)[:8]
    es = np.array(mc.get_excited(n=8, mode="ED"))
    assert np.max(np.abs(es - es_ref)) < ED_TOL
    assert es_ref[1] - es_ref[0] > 1e-3 # nondegenerate ground state
    e0, v0 = _gs(H)
    assert mc.gs_energy(mode="ED") == pytest.approx(e0, abs=ED_TOL)
    got = _observables(mc, mode="ED")
    for k, val in got.items():
        ref = np.vdot(v0, obs[k] @ v0)
        assert abs(val - ref) < 1e-9, k


@pytest.mark.parametrize("model", sorted(MODELS))
@pytest.mark.parametrize("v", BACKENDS)
def test_mixed_dmrg_matches_ed(model, v):
    labels, ham = MODELS[model]
    mc = _mc(labels, v)
    mc.nsweeps = 30
    mc.set_hamiltonian(ham(mc))
    e_ed = mc.gs_energy(mode="ED")
    obs_ed = _observables(mc, mode="ED")
    e_dmrg = mc.gs_energy(mode="DMRG")
    obs_dmrg = _observables(mc, mode="DMRG")
    assert e_dmrg == pytest.approx(e_ed, abs=E_TOL)
    for k in obs_ed:
        assert abs(obs_dmrg[k] - obs_ed[k]) < V_TOL, k


@pytest.mark.parametrize("v", BACKENDS)
def test_mixed_kpm_dynamical_correlator_dmrg_matches_ed(v):
    """S^zz(w) of the spin of the Kondo model, KPM on both solvers"""
    labels, ham = MODELS["kondo"]
    mc = _mc(labels, v)
    mc.set_hamiltonian(ham(mc))
    es = np.linspace(-0.5, 2.5, 120)
    kw = dict(name=(mc.Sz[1], mc.Sz[1]), es=es, delta=0.1, submode="KPM")
    x0, y0 = mc.get_dynamical_correlator(mode="ED", **kw)
    x1, y1 = mc.get_dynamical_correlator(mode="DMRG", **kw)
    assert len(y0) == len(y1) == len(es)
    scale = np.max(np.abs(y0.real))
    assert scale > 0.1 # a peak is in the window
    assert np.max(np.abs(y1 - y0))/scale < 2e-3


@pytest.mark.parametrize("v", BACKENDS)
def test_mixed_quench_dmrg_matches_ed(v):
    """<Sx_1(t)> of the Kondo spin after a quench: the ground state of the
    Kondo model with a transverse field on the spin, evolved without that
    field, TDVP on DMRG against expm_multiply on ED"""
    from dmrgpy import timedependent
    labels, ham = MODELS["kondo"]
    mc = _mc(labels, v)
    mc.tevol_method = "TDVP"
    out = {}
    for mode in ("ED", "DMRG"):
        mc.set_hamiltonian(ham(mc) + 0.3*mc.Sx[1])
        wf = mc.get_gs(mode=mode)
        mc.set_hamiltonian(ham(mc))
        _t, y = timedependent.evolve_and_measure(mc, operator=mc.Sx[1],
                nt=100, dt=0.1, wf=wf, mode=mode)
        out[mode] = np.asarray(y)
    assert len(out["ED"]) == len(out["DMRG"]) == 100
    assert np.ptp(out["ED"].real) > 0.05 # the spin really precesses
    assert np.max(np.abs(out["DMRG"] - out["ED"])) < QUENCH_TOL


# ---------------------------------------------------------------------
# Short chains: the automatic DMRG -> ED fallbacks
# ---------------------------------------------------------------------

def test_v3_short_chain_falls_back_to_ed():
    """One fermion location (ns=2) and two spins (ns=2) on v3"""
    mc = _mc(["F"], 3)
    assert mc.ns < 3
    mc.set_hamiltonian(1.3*mc.Nup[0]*mc.Ndn[0] - 0.4*mc.Ntot[0]
                       + 0.3*mc.Sx[0])
    # exact: empty 0, one electron -0.4 -+ 0.15, doubly occupied 0.5
    assert mc.gs_energy(mode="DMRG") == pytest.approx(-0.55, abs=ED_TOL)
    assert mc.vev(mc.Ntot[0]).real == pytest.approx(1.0, abs=ED_TOL)
    sc = _mc(["1/2", "1"], 3)
    sc.set_hamiltonian(sc.SS(0, 1))
    assert sc.gs_energy() == pytest.approx(-1.0, abs=ED_TOL) # S=1/2 doublet


def test_python_single_site_falls_back_to_ed():
    mc = _mc(["1"], "python")
    mc.set_hamiltonian(0.5*mc.Sz[0]*mc.Sz[0] - 0.2*mc.Sz[0])
    assert mc.gs_energy() == pytest.approx(0.0, abs=ED_TOL) # m=0
    mc.set_hamiltonian(0.5*mc.Sz[0]*mc.Sz[0] - 0.7*mc.Sz[0])
    assert mc.gs_energy() == pytest.approx(-0.2, abs=ED_TOL) # m=1


# ---------------------------------------------------------------------
# Conserved sectors
# ---------------------------------------------------------------------

def _sector_reference(H, obs, nf, sz2):
    """Lowest energy of the dense H restricted to the basis states with
    Nf=nf and 2Sz=sz2 (both charges are diagonal in the product basis)"""
    q1 = np.real(np.diag(obs["Nf"]))
    q2 = np.real(np.diag(obs["Sz2"]))
    sel = np.where((np.abs(q1-nf) < 1e-9) & (np.abs(q2-sz2) < 1e-9))[0]
    return np.linalg.eigvalsh(H[np.ix_(sel, sel)])[0]


@pytest.mark.parametrize("model,nf,sz2", [("kondo", 2, 1), ("kondo", 1, -2),
                                           ("spin1", 4, 0), ("spin1", 3, 1)])
def test_mixed_sector_ed_matches_restricted_oracle(model, nf, sz2):
    labels, ham = MODELS[model]
    mc = _mc(labels, mode="ED")
    mc.set_hamiltonian(ham(mc))
    mc.set_conserved_sector(Nf=nf, Sz=sz2)
    H, obs = _oracle(model)
    e_ref = _sector_reference(H, obs, nf, sz2)
    assert mc.gs_energy() == pytest.approx(e_ref, abs=ED_TOL)
    mc.set_conserved_sector(Sz=sz2) # Sz alone
    q2 = np.real(np.diag(obs["Sz2"]))
    sel = np.where(np.abs(q2-sz2) < 1e-9)[0]
    e_sz = np.linalg.eigvalsh(H[np.ix_(sel, sel)])[0]
    assert mc.gs_energy() == pytest.approx(e_sz, abs=ED_TOL)


@pytest.mark.parametrize("v", BACKENDS)
@pytest.mark.parametrize("sector", [dict(Nf=2), dict(Sz=1), dict(Nf=2, Sz=1)])
def test_mixed_sector_refused_on_dmrg(v, sector):
    """The DMRG sites of a mixed chain cannot carry total Sz or Nf. On an
    Ising-only Hamiltonian Nf+Sz used to be accepted by the session with
    Sz graded on the spin sites only."""
    mc = _mc(["F", "1/2", "F"], v)
    h = 0.5*mc.Sz[0]*mc.Sz[1] + 0.3*mc.Sz[1]*mc.Sz[2]
    hop = 0.4*(mc.Cdagup[0]*mc.Cup[2] + mc.Cdagdn[0]*mc.Cdn[2])
    mc.set_hamiltonian(h + hop + hop.get_dagger())
    with pytest.raises(NotImplementedError):
        mc.set_conserved_sector(**sector)
    assert mc.conserved_sector is None # the chain is left as it was
    assert np.isfinite(mc.gs_energy())


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("v", BACKENDS)
def test_mixed_sector_correlator_refused(v):
    """submode="SECTOR" sets Nf on an internal copy, which a chain with
    both site kinds cannot carry on DMRG: it is refused, and the message
    names the submodes that work"""
    mc = _mc(["F", "1/2", "F"], v)
    mc.set_hamiltonian(_kondo_mc(mc))
    with pytest.raises(NotImplementedError, match="submode=\"KPM\""):
        mc.get_dynamical_correlator(name=(mc.Sz[1], mc.Sz[1]),
                                    submode="SECTOR")


@pytest.mark.parametrize("v", BACKENDS)
def test_pure_sitesin_sectors_dmrg_match_ed(v):
    """All-spin sitesin conserves Sz on DMRG, all-"F" sitesin Nf; ED with
    the same sector gives the same energy, Sz on all-"F" is ED only"""
    sc = _mc(["1/2", "1", "1/2"], v)
    sc.set_hamiltonian(0.7*sc.SS(0, 1) + 1.3*sc.SS(1, 2) - 0.3*sc.Sz[1])
    sc.set_conserved_sector(Sz=2)
    assert sc.gs_energy(mode="DMRG") == pytest.approx(sc.gs_energy(mode="ED"),
                                                      abs=E_TOL)
    fc = _mc(["F", "F", "F"], v)
    h = 0
    for i in range(2):
        h = h + 0.8*(fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1])
    h = h + h.get_dagger()
    for i in range(3): h = h + 2.0*fc.Nup[i]*fc.Ndn[i] - 0.5*fc.Ntot[i]
    fc.set_hamiltonian(h)
    fc.set_conserved_sector(Nf=2)
    assert fc.gs_energy(mode="DMRG") == pytest.approx(fc.gs_energy(mode="ED"),
                                                      abs=E_TOL)
    with pytest.raises(NotImplementedError):
        fc.set_conserved_sector(Sz=0)
