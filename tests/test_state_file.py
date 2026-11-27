"""Wavefunctions on disk (statefile.py): wf.save(name),
chain.save_state(wf, name) and chain.load_state(name), on
itensor_version=3, "python", "julia_live" and ED."""
import itertools
import os
import subprocess
import sys

import numpy as np
import pytest

from dmrgpy import spinchain, fermionchain, bosonchain, parafermionchain
from dmrgpy import mixedchain
from dmrgpy import cppext, statefile
from _helpers import julia_live_param

V3 = pytest.param(3, marks=pytest.mark.skipif(not cppext.available(3),
        reason="needs the compiled mpscpp3 extension"))
BACKENDS = ["python", V3, julia_live_param()]
SESSION_BACKENDS = ["python", V3]
SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(
        __file__))), "src")


# ----------------------------------------------------------------- chains
# Hamiltonians with no symmetry to hide a basis reordering, a sign of a
# fermion basis state or a conjugation: fields along x, y and z, complex
# hoppings, pairing.

def _settle(c, h):
    c.maxm, c.nsweeps = 30, 10
    c.set_hamiltonian(h)


def spin_chain(iv, sites=(2,)*6):
    sc = spinchain.Spin_Chain(list(sites), itensor_version=iv)
    r = np.random.default_rng(1)
    h, obs = 0, []
    for i in range(sc.ns):
        h = h + r.normal()*sc.Sx[i] + r.normal()*sc.Sy[i] + r.normal()*sc.Sz[i]
        obs += [sc.Sx[i], sc.Sy[i], sc.Sz[i]]
    for i in range(sc.ns-1):
        h = (h + sc.Sx[i]*sc.Sx[i+1] + r.normal()*sc.Sy[i]*sc.Sz[i+1]
             + sc.Sz[i]*sc.Sz[i+1])
        obs += [sc.Sx[i]*sc.Sy[i+1]]
    _settle(sc, h)
    return sc, h, obs


def mixed_spin_chain(iv):
    return spin_chain(iv, sites=(3, 4, 5, 6, 2))


def _fermion_model(fc, C, Cd):
    r = np.random.default_rng(2)
    m = len(C)
    h, obs = 0, []
    for i in range(m):
        h = h + r.normal()*Cd[i]*C[i]
        for j in range(i+1, m):
            t = r.normal() + 1j*r.normal()
            d = 0.3*(r.normal() + 1j*r.normal())
            h = h + t*Cd[i]*C[j] + np.conj(t)*Cd[j]*C[i]
            h = h + d*Cd[i]*Cd[j] + np.conj(d)*C[j]*C[i]
            obs += [Cd[i]*C[j], C[i]*C[j]]
    h = h + 2.0*Cd[0]*C[0]*Cd[1]*C[1]
    _settle(fc, h)
    return fc, h, obs


def spinless_chain(iv):
    fc = fermionchain.Fermionic_Chain(5, itensor_version=iv)
    return _fermion_model(fc, list(fc.C), list(fc.Cdag))


def native_spinful_chain(iv):
    fc = fermionchain.Spinful_Fermionic_Chain_Native(3, itensor_version=iv)
    C = [c for i in range(3) for c in (fc.Cup[i], fc.Cdn[i])]
    Cd = [c for i in range(3) for c in (fc.Cdagup[i], fc.Cdagdn[i])]
    return _fermion_model(fc, C, Cd)


def boson_chain(iv):
    bc = bosonchain.Bosonic_Chain(4, maxnb=[3]*4, itensor_version=iv)
    r = np.random.default_rng(4)
    h, obs = 0, []
    for i in range(4):
        h = (h + r.normal()*bc.N[i] + 0.5*bc.N[i]*bc.N[i]
             + 0.3*(bc.A[i] + bc.Adag[i]))
        obs += [bc.N[i], bc.A[i]]
    for i in range(3):
        t = r.normal() + 1j*r.normal()
        h = h + t*bc.Adag[i]*bc.A[i+1] + np.conj(t)*bc.Adag[i+1]*bc.A[i]
        obs += [bc.Adag[i]*bc.A[i+1]]
    _settle(bc, h)
    return bc, h, obs


def parafermion_chain(iv):
    pf = parafermionchain.Parafermionic_Chain(4, itensor_version=iv)
    r = np.random.default_rng(5)
    h, obs = 0, []
    for i in range(3):
        z = r.normal() + 1j*r.normal()
        h = h + z*pf.Chi[i]*pf.Chid[i+1]
        obs += [pf.Chi[i]*pf.Chid[i+1]]
    h = h + h.get_dagger()
    _settle(pf, h)
    return pf, h, obs


def mixed_spin_fermion_chain(iv):
    """A fermion location between two spins: the hopping string of the
    two fermion locations crosses a spin site."""
    mc = mixedchain.Mixed_Spin_Fermion_Chain(["F", "1/2", "F", "1"],
                                             itensor_version=iv)
    r = np.random.default_rng(6)
    h, obs = 0, []
    for i in range(4):
        h = h + r.normal()*mc.Sx[i] + r.normal()*mc.Sy[i] + r.normal()*mc.Sz[i]
        obs += [mc.Sx[i], mc.Sy[i], mc.Sz[i]]
    for i in range(3):
        h = h + mc.Sx[i]*mc.Sx[i+1] + mc.Sz[i]*mc.Sz[i+1]
    t = r.normal() + 1j*r.normal()
    hop = t*(mc.Cdagup[0]*mc.Cup[2] + 0.5*mc.Cdagdn[0]*mc.Cdn[2])
    h = h + hop + hop.get_dagger() + 0.7*mc.Ntot[0] - 0.4*mc.Ntot[2]
    obs += [mc.Cdagup[0]*mc.Cup[2], mc.Cdagdn[0]*mc.Cdn[2], mc.Ntot[0],
            mc.Ntot[2]]
    _settle(mc, h)
    return mc, h, obs


CHAINS = {"spin": spin_chain, "mixed_spin_fermion": mixed_spin_fermion_chain, "mixed_spins": mixed_spin_chain,
          "spinless": spinless_chain, "native_spinful": native_spinful_chain,
          "boson": boson_chain, "parafermion": parafermion_chain}


def _heisenberg(iv, n=6, **kwargs):
    sc = spinchain.Spin_Chain([2]*n, itensor_version=iv, **kwargs)
    h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
            for i in range(n-1))
    _settle(sc, h)
    return sc, h


def _values(c, ops, wf):
    return np.array([c.vev(o, wf=wf) for o in ops])


def _bond_dimensions(c, wf):
    if hasattr(wf, "jlmps"):
        from dmrgpy.mpsjulialive.juliasession import Main
        arrs = Main.mps_site_arrays(wf.jlmps)
    else:
        arrs = c._session.mps_site_arrays(wf.cpp_handle)
    return [np.asarray(a).shape[2] for a in arrs]


# ---------------------------------------------------------- the round trip

@pytest.mark.parametrize("iv", BACKENDS)
def test_round_trip_keeps_the_state_and_its_norm(iv, tmp_path):
    """A complex state that is not normalized comes back as it was: the
    same norm, and an overlap with the original equal to the norm squared,
    which a conjugation of the arrays would turn into another number."""
    sc, h, obs = spin_chain(iv)
    wf = (0.7-0.4j)*sc.get_gs()
    fname = wf.save(tmp_path/"state")
    assert fname == str(tmp_path/"state.npz") and os.path.isfile(fname)
    assert os.listdir(tmp_path) == ["state.npz"]  # no temporary file left
    back = sc.load_state(tmp_path/"state")
    n2 = wf.dot(wf)
    assert n2.real == pytest.approx(0.65, abs=1e-10)
    assert back.dot(back) == pytest.approx(n2, abs=1e-10)
    assert wf.dot(back) == pytest.approx(n2, abs=1e-10)
    assert np.max(np.abs(_values(sc, [h]+obs, wf)
                         - _values(sc, [h]+obs, back))) < 1e-10
    assert _bond_dimensions(sc, back) == _bond_dimensions(sc, wf)
    assert max(_bond_dimensions(sc, wf)) > 1
    # the two spellings of the name are one file, and the chain-level
    # writer writes the same thing
    again = sc.load_state(str(tmp_path/"state.npz"))
    assert wf.dot(again) == pytest.approx(n2, abs=1e-10)
    sc.save_state(wf, tmp_path/"state2.npz")
    assert wf.dot(sc.load_state(tmp_path/"state2")) == pytest.approx(
            n2, abs=1e-10)


@pytest.mark.parametrize("model", sorted(CHAINS))
def test_file_of_one_backend_loads_on_the_other(model, tmp_path):
    """The site arrays are in the basis of the site types, which v3 and
    "python" share: every observable of a low-symmetry state agrees after
    the state crosses, in both directions."""
    if not cppext.available(3):
        pytest.skip("needs the compiled mpscpp3 extension")
    chains = {iv: CHAINS[model](iv) for iv in ("python", 3)}
    for a, b in itertools.permutations(chains, 2):
        (ca, ha, oa), (cb, hb, ob) = chains[a], chains[b]
        wf = (0.7-0.4j)*ca.get_gs()
        wf.save(tmp_path/"x")
        assert statefile.read_header(tmp_path/"x")["backend"] == str(a)
        back = cb.load_state(tmp_path/"x")
        assert back.dot(back) == pytest.approx(wf.dot(wf), abs=1e-10)
        assert np.max(np.abs(_values(ca, [ha]+oa, wf)
                             - _values(cb, [hb]+ob, back))) < 1e-9


@pytest.mark.parametrize("iv", [julia_live_param()])
@pytest.mark.parametrize("model", ["spin", "spinless"])
def test_file_crosses_between_julia_and_python(iv, model, tmp_path):
    chains = {v: CHAINS[model](v) for v in ("python", iv)}
    for a, b in itertools.permutations(chains, 2):
        (ca, ha, oa), (cb, hb, ob) = chains[a], chains[b]
        wf = (0.7-0.4j)*ca.get_gs()
        wf.save(tmp_path/"x")
        back = cb.load_state(tmp_path/"x")
        assert back.dot(back) == pytest.approx(wf.dot(wf), abs=1e-10)
        assert np.max(np.abs(_values(ca, [ha]+oa, wf)
                             - _values(cb, [hb]+ob, back))) < 1e-9


# ------------------------------------------ what the file is for: a new job

_WRITER = """
import sys
sys.path.insert(0, %r)
import numpy as np
from dmrgpy import spinchain
sc = spinchain.Spin_Chain([2]*8, itensor_version=%r)
sc.maxm, sc.nsweeps = 30, 10
sc.set_hamiltonian(sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
        + sc.Sz[i]*sc.Sz[i+1] + 0.3*(-1)**i*sc.Sz[i] for i in range(7)))
e = sc.gs_energy()
sc.get_gs().save(%r)
print(repr(float(e)), repr(float(sc.vev(sc.Sz[2]).real)))
"""


@pytest.mark.parametrize("iv", SESSION_BACKENDS)
def test_state_written_by_another_process(iv, tmp_path):
    """The use case: one process solves and saves, a later one builds the
    same chain, loads the state and makes it the chain's. The next
    gs_energy() is <H> of that state, and a vev measures it."""
    name = str(tmp_path/"gs")
    out = subprocess.run([sys.executable, "-c", _WRITER % (SRC, iv, name)],
                         capture_output=True, text=True, check=True)
    e_ref, sz_ref = (float(x) for x in out.stdout.split()[-2:])
    sc = spinchain.Spin_Chain([2]*8, itensor_version=iv)
    sc.maxm, sc.nsweeps = 30, 10
    sc.set_hamiltonian(sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
            + sc.Sz[i]*sc.Sz[i+1] + 0.3*(-1)**i*sc.Sz[i] for i in range(7)))
    sc.set_gs(sc.load_state(name))
    assert sc.gs_energy() == pytest.approx(e_ref, abs=1e-9)
    assert sc.vev(sc.Sz[2]).real == pytest.approx(sz_ref, abs=1e-9)


@pytest.mark.parametrize("iv", BACKENDS)
def test_set_gs_of_a_loaded_state_is_not_swept(iv, tmp_path):
    """A loaded state that is not the ground state stays what it is after
    set_gs(): the energy is its <H>, above E0, and vev reads it."""
    sc, h = _heisenberg(iv)
    e0 = sc.gs_energy()
    wf = sc.applyoperator(sc.Sz[1], sc.get_gs())
    wf = wf*(1./wf.norm())
    e_wf = sc.vev(h, wf=wf).real
    assert e_wf > e0 + 0.1
    wf.save(tmp_path/"excited")
    fresh, _ = _heisenberg(iv)
    fresh.set_gs(fresh.load_state(tmp_path/"excited"))
    assert fresh.gs_energy() == pytest.approx(e_wf, abs=1e-8)
    assert fresh.vev(h).real == pytest.approx(e_wf, abs=1e-8)


# ----------------------------------------------------------------- sectors

@pytest.mark.parametrize("iv", BACKENDS)
def test_sector_state_is_written_as_its_dense_form(iv, tmp_path):
    """A state solved in a sector is stored as the dense state
    promote_mps() gives, and loads on the chain with the sector cleared."""
    sc, h = _heisenberg(iv)
    sc.set_conserved_sector(Sz=2)  # in units of 1/2: total Sz = 1
    e = sc.gs_energy()
    wf = sc.get_gs()
    wf.save(tmp_path/"sz1")
    assert statefile.read_header(tmp_path/"sz1")["sector"] == {"Sz": 2}
    with pytest.raises(ValueError, match="dense sites only"):
        sc.load_state(tmp_path/"sz1")
    dense, _ = _heisenberg(iv)
    back = dense.load_state(tmp_path/"sz1")
    assert back.dot(back).real == pytest.approx(1.0, abs=1e-10)
    assert dense.vev(h, wf=back).real == pytest.approx(e, abs=1e-8)
    sztot = sum(dense.Sz[i] for i in range(dense.ns))
    assert dense.vev(sztot, wf=back).real == pytest.approx(1.0, abs=1e-8)
    promoted = sc.promote_mps(wf)
    sc.set_conserved_sector()
    assert abs(sc.overlap(promoted, sc.load_state(tmp_path/"sz1"))) \
        == pytest.approx(1.0, abs=1e-10)


# ---------------------------------------------------------------------- ED

def test_ed_state_round_trip(tmp_path):
    sc, h = _heisenberg("python")
    wf = (0.5+0.2j)*sc.get_gs(mode="ED")
    sc.save_state(wf, tmp_path/"ed")
    meta = statefile.read_header(tmp_path/"ed")
    assert meta["kind"] == "ED" and meta["dim"] == 64
    assert meta["sites"] == ["2"]*6
    fresh, _ = _heisenberg("python")
    back = fresh.load_state(tmp_path/"ed")
    assert back.mode == "ED"
    assert np.allclose(back.v, wf.v, atol=0, rtol=0)
    # a state that is not the ground state, so that a chain which ignored
    # the set state would answer with another number
    up = sc.Sz[1]*sc.get_gs(mode="ED")
    sc.save_state(up, tmp_path/"ed_up")
    back = fresh.load_state(tmp_path/"ed_up")
    e_up = (np.vdot(up.v, (h*up).v)/np.vdot(up.v, up.v)).real
    assert e_up > sc.gs_energy(mode="ED") + 0.1
    assert fresh.vev(h, wf=back, mode="ED").real/np.vdot(up.v, up.v).real \
        == pytest.approx(e_up, abs=1e-10)
    # the state's own writer has no site list to record
    wf.save(tmp_path/"ed_bare")
    assert statefile.read_header(tmp_path/"ed_bare")["sites"] is None
    assert np.allclose(fresh.load_state(tmp_path/"ed_bare").v, wf.v)


def test_ed_state_of_a_sector_needs_the_same_sector(tmp_path):
    sc, h = _heisenberg("python")
    sc.set_conserved_sector(Sz=0)
    wf = sc.get_gs(mode="ED")
    assert wf.v.size == 20
    sc.save_state(wf, tmp_path/"ed0")
    assert np.allclose(sc.load_state(tmp_path/"ed0").v, wf.v)
    other, _ = _heisenberg("python")
    with pytest.raises(ValueError, match="sector"):
        other.load_state(tmp_path/"ed0")
    other.set_conserved_sector(Sz=2)
    with pytest.raises(ValueError, match="sector"):
        other.load_state(tmp_path/"ed0")


# ---------------------------------------------------------------- refusals

def test_a_file_of_another_chain_is_refused(tmp_path):
    sc, h = _heisenberg("python")
    sc.get_gs().save(tmp_path/"six")
    longer, _ = _heisenberg("python", n=8)
    with pytest.raises(ValueError, match="sites"):
        longer.load_state(tmp_path/"six")
    spin1 = spinchain.Spin_Chain([3]*6, itensor_version="python")
    with pytest.raises(ValueError, match="sites"):
        spin1.load_state(tmp_path/"six")
    fc = fermionchain.Fermionic_Chain(6, itensor_version="python")
    with pytest.raises(ValueError, match="sites"):
        fc.load_state(tmp_path/"six")


def test_ed_file_on_a_fermion_chain_with_no_hamiltonian(tmp_path):
    """The ED object of a fermion chain is built from its Hamiltonian;
    the refusal names set_hamiltonian()."""
    fc = fermionchain.Fermionic_Chain(4, itensor_version="python")
    fc.set_hamiltonian(sum(fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
                           for i in range(3)))
    wf = fc.get_gs(mode="ED")
    fc.save_state(wf, tmp_path/"fed")
    assert np.allclose(fc.load_state(tmp_path/"fed").v, wf.v)
    bare = fermionchain.Fermionic_Chain(4, itensor_version="python")
    with pytest.raises(RuntimeError, match="set_hamiltonian"):
        bare.load_state(tmp_path/"fed")


def test_a_file_that_is_not_a_state_is_refused(tmp_path):
    sc, h = _heisenberg("python")
    np.savez(tmp_path/"other.npz", x=np.arange(3))
    with pytest.raises(ValueError, match="not a dmrgpy state file"):
        sc.load_state(tmp_path/"other")
    with pytest.raises(FileNotFoundError):
        sc.load_state(tmp_path/"missing")


def test_a_newer_format_is_refused(tmp_path, monkeypatch):
    sc, h = _heisenberg("python")
    monkeypatch.setattr(statefile, "VERSION", statefile.VERSION+1)
    sc.get_gs().save(tmp_path/"new")
    monkeypatch.undo()
    with pytest.raises(ValueError, match="format version"):
        sc.load_state(tmp_path/"new")


def test_v2_has_no_array_bridge(tmp_path):
    if not cppext.available(2):
        pytest.skip("needs the compiled mpscpp2 extension")
    sc, h = _heisenberg(2)
    with pytest.raises(NotImplementedError, match="itensor_version"):
        sc.get_gs().save(tmp_path/"v2")
    py, _ = _heisenberg("python")
    py.get_gs().save(tmp_path/"py")
    with pytest.raises(NotImplementedError, match="itensor_version"):
        sc.load_state(tmp_path/"py")
    assert not os.path.exists(tmp_path/"v2.npz")
