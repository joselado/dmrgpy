"""The speedups of the ED Kondo route (docs/kondo_speedup_plan.md, items 1
to 6). Every test pins the numbers against the form of the code before the
change, kept here, and the saving, as a count of F or Theta entries, of
coefficient entries or of diagonalizations, so that undoing the saving fails
a test."""
import numpy as np
import pytest

from dmrgpy import spinchain
from dmrgpy.algebra import algebra
from dmrgpy.kondospectrumtk import conductance as cd
from dmrgpy.kondospectrumtk.edkondo import KondoSpectrum

_EPS3 = np.zeros((3, 3, 3))
for _i, _j, _k in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
    _EPS3[_i, _j, _k] = 1.
for _i, _j, _k in [(0, 2, 1), (2, 1, 0), (1, 0, 2)]:
    _EPS3[_i, _j, _k] = -1.


def _old_coefficients(ks, occ):
    """Im of the Levi-Civita triple product over occ x every f x every m,
    in complex arithmetic, as built before the speedups."""
    Xi = np.stack([np.asarray(ks.Sx, dtype=complex), np.asarray(ks.Sy, dtype=complex),
                   np.asarray(ks.Sz, dtype=complex)], axis=-1)
    return np.imag(np.einsum('jkl,ifl,fmk,mij->ifm', _EPS3, Xi[occ], Xi,
                             Xi[:, occ], optimize=True))


def _old_second(ks, eVs, T0=1.0, U=0.0):
    """Second order over the occupied initial states and every final state."""
    eVs = np.asarray(eVs, dtype=float)
    occ = cd._occupied_states(ks)
    S = [np.asarray(a, dtype=complex)[:, occ] for a in (ks.Sx, ks.Sy, ks.Sz)]
    M2 = (0.5*np.abs(S[0] - 1j*S[1])**2 + 0.5*np.abs(S[0] + 1j*S[1])**2
          + np.abs(S[2])**2)
    eps = ks.e[:, None] - ks.e[occ][None, :]
    w = ks.p[occ][None, :]*M2
    spin = (np.einsum('fi,efi->e', w, cd._theta_raw(ks, eVs[:, None, None] - eps))
            + np.einsum('fi,efi->e', w, cd._theta_raw(ks, -eVs[:, None, None] - eps)))
    return 2*np.pi*T0**2*(spin + 4*U**2)


def _old_third(ks, eVs, Jrho_s, T0=1.0, omega0=20e-3, Gamma0=5e-6, Fb=None):
    """Third-order Kondo term in one pass over every f and m, both F tables
    of both diagrams evaluated in full, as it was before the speedups."""
    eVs = np.asarray(eVs, dtype=float)
    occ = cd._occupied_states(ks)
    p = ks.p[occ]
    coeff = _old_coefficients(ks, occ)/2.
    eps = ks.e[None, :] - ks.e[occ, None]
    e_fm = ks.e[:, None] - ks.e[None, :]
    Fcall = cd._get_F(ks, omega0, Gamma0, Fb)
    dim, nocc = ks.dim, len(occ)

    def one(v):
        nv = len(v)
        Th = cd._theta_raw(ks, v[:, None, None] - eps[None])
        Fd = Fcall((v[:, None, None] - eps[None]).ravel()).reshape(nv, nocc, dim)
        Fx = Fcall((v[:, None, None] - e_fm[None]).ravel()).reshape(nv, dim, dim)
        return (np.einsum('i,ifm,eif,eim->e', p, coeff, Th, Fd, optimize=True)
                + np.einsum('i,ifm,eif,efm->e', p, coeff, Th, Fx, optimize=True))
    return 4*np.pi*T0**2*Jrho_s*(one(eVs) + one(-eVs))


def _old_potential(ks, eVs, Jrho_s, U, T0=1.0, omega0=20e-3, Gamma0=5e-6,
                   Fb=None):
    """Third-order potential term with four F arrays and no bias chunks."""
    eVs = np.asarray(eVs, dtype=float)
    occ = cd._occupied_states(ks)
    p = ks.p[occ]
    Xi = np.stack([np.asarray(a, dtype=complex) for a in (ks.Sx, ks.Sy, ks.Sz)], axis=-1)
    loop = np.real(np.einsum('imk,mik->im', Xi[occ], Xi[:, occ]))
    eps = ks.e[None, :] - ks.e[occ, None]
    Fcall = cd._get_F(ks, omega0, Gamma0, Fb)
    nocc, dim = len(occ), ks.dim

    def one(v):
        Fim = Fcall((v[:, None, None] - eps[None]).ravel()).reshape(len(v), nocc, dim)
        Fmi = Fcall((v[:, None, None] + eps[None]).ravel()).reshape(len(v), nocc, dim)
        return np.einsum('i,im,eim->e', p, loop, Fim - Fmi)*cd._theta_raw(ks, v)
    return 4*np.pi*T0**2*Jrho_s*U*(one(eVs) - one(-eVs))


def _chain(n=5, field=True, complex_h=False):
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    h = 0
    for i in range(n - 1):
        h = h + 3e-3*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
                      + 0.7*sc.Sz[i]*sc.Sz[i+1])
    if field: h = h + 0.4e-3*sc.Sz[0] + 0.3e-3*sc.Sx[2]
    if complex_h: h = h + 0.5e-3*sc.Sy[1]
    sc.set_hamiltonian(h)
    return sc


EV = np.linspace(-2e-3, 2e-3, 9)


def _rel(a, b):
    return np.max(np.abs(np.asarray(a) - np.asarray(b)))/np.max(np.abs(b))


class _CountF:
    """A callable F that counts the entries it is asked for."""
    def __init__(self, f):
        self.f, self.n = f, 0
        for name in ("odd_part",):
            if hasattr(f, name): setattr(self, name, getattr(f, name))
    def __call__(self, x):
        self.n += np.size(x)
        return self.f(x)


cd_F0 = cd.F0


def _patch_F0(monkeypatch):
    calls = []
    def counting(x, **kw):
        calls.append(np.size(x))
        return cd_F0(x, **kw)
    monkeypatch.setattr(cd, "F0", counting)
    return calls


# item 1: final states closed by the step


@pytest.mark.parametrize("T", [0.0, 0.5])
def test_final_states_cut_keeps_the_numbers(T):
    ks = KondoSpectrum(_chain(), 0, T)
    occ = cd._occupied_states(ks)
    fin = cd._final_states(ks, occ, EV)
    assert len(fin) < ks.dim # the cut is not empty on this chain
    tol = 1e-13 if T == 0. else 1e-12
    assert _rel(cd.second_order_dIdV(ks, EV, U=0.2), _old_second(ks, EV, U=0.2)) < tol
    Fb = cd._get_F(ks, 20e-3, 5e-6, None)
    assert _rel(cd.third_order_kondo_dIdV(ks, EV, -0.05, Fb=Fb if T else None),
                _old_third(ks, EV, -0.05, Fb=Fb if T else None)) < tol


def test_final_states_cut_at_zero_temperature_drops_only_closed_steps():
    """Every state left out has e_f - e_i > max|eV| for every occupied i, so
    its Theta0 is exactly 0; and the T=0 F count is that of the kept ones."""
    ks = KondoSpectrum(_chain(), 0, 0.)
    occ = cd._occupied_states(ks)
    fin = cd._final_states(ks, occ, EV)
    out = np.setdiff1d(np.arange(ks.dim), fin)
    eps = ks.e[None, out] - ks.e[occ, None]
    assert np.all(cd._theta_raw(ks, np.max(np.abs(EV)) - eps) == 0.)
    assert np.all(cd._theta_raw(ks, np.max(np.abs(EV)) - (ks.e[None, fin]
                  - ks.e[occ, None])).max(axis=0) > 0.)


def test_final_states_cut_counts_F_entries(monkeypatch):
    """The saving: n_f*dim exchange entries per bias and direction, not dim^2."""
    ks = KondoSpectrum(_chain(), 0, 0.)
    occ = cd._occupied_states(ks)
    nf, nocc, dim = len(cd._final_states(ks, occ, EV)), len(occ), ks.dim
    calls = _patch_F0(monkeypatch)
    cd.third_order_kondo_dIdV(ks, EV, -0.05)
    assert sum(calls) <= 2*len(EV)*(nocc + nf)*dim
    assert sum(calls) < 2*len(EV)*(nocc + dim)*dim


def test_final_states_cut_on_an_unsorted_spectrum():
    """The TASK and Krylov spectra need not be sorted: the cut is a mask."""
    ks = KondoSpectrum(_chain(), 0, 0.5)
    ref2 = cd.second_order_dIdV(ks, EV)
    ref3 = cd.third_order_kondo_dIdV(ks, EV, -0.05)
    perm = np.random.default_rng(3).permutation(ks.dim)
    ks.e, ks.p = ks.e[perm], ks.p[perm]
    for name in ("Sx", "Sy", "Sz"):
        setattr(ks, name, getattr(ks, name)[np.ix_(perm, perm)])
    assert _rel(cd.second_order_dIdV(ks, EV), ref2) < 1e-13
    assert _rel(cd.third_order_kondo_dIdV(ks, EV, -0.05), ref3) < 1e-12


# item 3: real arithmetic for a real H


class _OldSpectrum:
    """KondoSpectrum as it was built before item 3: the complex eigenpairs of
    the ED object and complex transforms of the spin operators."""
    def __init__(self, chain, site, T, kB=8.617333262e-5, degeneracy_tol=1e-9):
        edobj = chain.get_ED_obj()
        emu, vs = edobj.get_diagonalized_hamiltonian()
        emu = np.array(emu, dtype=float)
        order = np.argsort(emu)
        emu, vs = emu[order], np.array(vs)[:, order]
        self.e, self.dim, self.T, self.kB = emu - emu[0], len(emu), T, kB
        if T == 0.:
            gs = (self.e <= degeneracy_tol*max(self.e[-1], 1e-300)).astype(float)
            self.p = gs/np.sum(gs)
        else:
            p = np.exp(-self.e/(kB*T))
            self.p = p/np.sum(p)
        for name in ("Sx", "Sy", "Sz"):
            op = getattr(chain, name)[site]
            mat = np.array(algebra.todense(edobj.get_operator(op)), dtype=complex)
            setattr(self, name, np.conjugate(vs.T)@mat@vs)


def _all_terms(ks, Fb=None):
    return (cd.second_order_dIdV(ks, EV, U=0.2),
            cd.third_order_kondo_dIdV(ks, EV, -0.05, Fb=Fb),
            cd.third_order_potential_dIdV(ks, EV, -0.05, 0.2, Fb=Fb))


def test_real_hamiltonian_gives_real_arrays():
    """The saving: real eigenvectors, real Sx and Sz, Sy = i times a real
    matrix, and the real form read by the Levi-Civita build. A complex H
    (a field along y) keeps the complex path."""
    ks = KondoSpectrum(_chain(), 0, 0.5)
    assert ks.vs.dtype == np.float64
    assert ks.Sx.dtype == np.float64 and ks.Sz.dtype == np.float64
    assert np.iscomplexobj(ks.Sy) and not np.any(ks.Sy.real)
    assert cd._spin_triplet(ks)[1]
    kc = KondoSpectrum(_chain(complex_h=True), 0, 0.5)
    assert np.iscomplexobj(kc.vs) and not cd._spin_triplet(kc)[1]


def test_real_levi_civita_equals_the_complex_product():
    ks = KondoSpectrum(_chain(), 0, 0.5)
    occ = cd._occupied_states(ks)
    fin = cd._final_states(ks, occ, EV)
    ref = np.imag(cd._triple_product_coefficients(ks, occ, fin))
    S, real = cd._spin_triplet(ks)
    got = cd._levi_civita(S, real, occ, fin)
    assert got.dtype == np.float64
    assert _rel(got, ref) < 1e-13
    kc = KondoSpectrum(_chain(complex_h=True), 0, 0.5)
    ref = np.imag(cd._triple_product_coefficients(kc, occ, fin))
    assert _rel(cd._levi_civita(*cd._spin_triplet(kc), occ, fin), ref) < 1e-13


@pytest.mark.parametrize("complex_h", [False, True])
def test_numba_and_numpy_builds_agree(complex_h):
    """The kernel above NUMBA_ENTRIES and NumPy below it make the same
    operations in the same order: the same bits on the real form."""
    ks = KondoSpectrum(_chain(complex_h=complex_h), 0, 2.)
    rows, cols = np.arange(0, ks.dim, 3), np.arange(ks.dim)[::-1]
    S, real = cd._spin_triplet(ks)
    a = cd._levi_civita(S, real, rows, cols, numba=True)
    b = cd._levi_civita(S, real, rows, cols, numba=False)
    if real: assert np.array_equal(a, b)
    else: assert _rel(a, b) < 1e-14


@pytest.mark.parametrize("field", [True, False])
@pytest.mark.parametrize("T", [0.0, 2.0])
def test_real_path_equals_complex_path(field, T):
    """Without a field the 5-spin chain has a doublet ground state and
    degenerate multiplets, where the real and complex eigh give different
    bases; every sum is a trace over each level and does not see it."""
    new = KondoSpectrum(_chain(field=field), 0, T)
    old = _OldSpectrum(_chain(field=field), 0, T)
    if not field: assert np.sum(new.p > 0.4*np.max(new.p)) >= 2
    Fb = cd._get_F(new, 20e-3, 5e-6, None) if T else None
    for a, b in zip(_all_terms(new, Fb), _all_terms(old, Fb)):
        assert _rel(a, b) < 1e-12


def test_complex_hamiltonian_keeps_the_numbers():
    new = KondoSpectrum(_chain(complex_h=True), 0, 2.0)
    old = _OldSpectrum(_chain(complex_h=True), 0, 2.0)
    Fb = cd._get_F(new, 20e-3, 5e-6, None)
    for a, b in zip(_all_terms(new, Fb), _all_terms(old, Fb)):
        assert _rel(a, b) < 1e-12


def _impurity_chain(complex_h):
    """An S=1 impurity on three spins 1/2: its sparse H is not in canonical
    form, which is where reading H.imag in place once scrambled H."""
    sc = spinchain.Spin_Chain(["S=1"] + ["S=1/2"]*3)
    h = 1e-3*sc.Sz[0]*sc.Sz[0] + 2e-4*sc.Sz[0]
    for i in range(3):
        h = h + 1e-3*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
    if complex_h: h = h + 3e-4*sc.Sy[0] + 1e-4*sc.Sx[1]
    sc.set_hamiltonian(h)
    return sc


@pytest.mark.parametrize("complex_h", [False, True])
def test_the_hamiltonian_matrix_is_not_touched(complex_h):
    sc = _impurity_chain(complex_h)
    H = sc.get_ED_obj().get_hamiltonian()
    before = (H.data.copy(), H.indices.copy(), H.indptr.copy())
    new = KondoSpectrum(sc, 0, 0.5)
    for a, b in zip(before, (H.data, H.indices, H.indptr)):
        assert np.array_equal(a, b)
    old = _OldSpectrum(_impurity_chain(complex_h), 0, 0.5)
    assert np.max(np.abs(new.e - old.e)) < 1e-15
    for a, b in zip(_all_terms(new, cd._get_F(new, 20e-3, 5e-6, None)),
                    _all_terms(old, cd._get_F(old, 20e-3, 5e-6, None))):
        assert _rel(a, b) < 1e-12


def test_real_eigenpairs_are_kept_and_the_ed_object_is_untouched(monkeypatch):
    """Two spectra of one chain diagonalize once, and the ED object's own
    (complex) diagonalization, which edtk/dynamics.py reads, is not changed."""
    calls = []
    eigh = algebra.eigh
    monkeypatch.setattr(algebra, "eigh", lambda m: (calls.append(1), eigh(m))[1])
    sc = _chain()
    a = KondoSpectrum(sc, 0, 0.)
    b = KondoSpectrum(sc, 0, 2.)
    assert len(calls) == 1
    assert np.array_equal(a.vs, b.vs)
    edobj = sc.get_ED_obj()
    assert edobj.Diagonalized_Hamiltonian is None
    emu, vs = edobj.get_diagonalized_hamiltonian()
    assert np.iscomplexobj(vs) and len(calls) == 2
    sc.set_hamiltonian(sc.hamiltonian + 1e-4*sc.Sz[1]) # a new H: a new eigh
    KondoSpectrum(sc, 0, 0.)
    assert len(calls) == 3


# item 5: the potential term in chunks of biases, two F arrays by evenness


def test_fbuilder_is_even_to_roundoff():
    """What the even forms rest on at T>0: the table of FBuilder is even to
    its roundoff, which it measures at construction (odd_part)."""
    for T in (0.5, 2.0, 20.0):
        fb = cd.FBuilder(T, omega0=20e-3, Gamma0=5e-6)
        assert fb.odd_part < cd.EVEN_TOL
    assert cd.FBuilder(2.0, tabulate=False).odd_part is None


def test_potential_term_halves_F_at_zero_temperature(monkeypatch):
    ks = KondoSpectrum(_chain(), 0, 0.)
    nocc = len(cd._occupied_states(ks))
    ref = _old_potential(ks, EV, -0.05, 0.2)
    calls = _patch_F0(monkeypatch)
    got = cd.third_order_potential_dIdV(ks, EV, -0.05, 0.2)
    assert _rel(got, ref) < 1e-14
    assert sum(calls) == 2*len(EV)*nocc*ks.dim # four arrays before


@pytest.mark.parametrize("trust", [True, False])
def test_potential_term_at_finite_temperature(trust):
    """The even form only when the FBuilder says its table is even; a
    callable without odd_part keeps the four arrays."""
    ks = KondoSpectrum(_chain(), 0, 2.0)
    nocc = len(cd._occupied_states(ks))
    fb = cd.FBuilder(2.0, omega0=20e-3, Gamma0=5e-6)
    ref = _old_potential(ks, EV, -0.05, 0.2, Fb=fb)
    F = _CountF(fb)
    if not trust: del F.odd_part
    got = cd.third_order_potential_dIdV(ks, EV, -0.05, 0.2, Fb=F)
    assert _rel(got, ref) < (1e-12 if trust else 1e-14)
    assert F.n == (2 if trust else 4)*len(EV)*nocc*ks.dim


def test_potential_term_memory_is_chunked(monkeypatch):
    """Every state occupied (8 spins at 1000 K) and 61 biases: the old term
    made its four (61, 256, 256) arrays at once, about 130 MB. The exact
    form is used here: at 1000 K the term is a near cancellation, 1e-7 of
    its scale, and the 1e-14 odd part of the table is 2e-11 of the term."""
    import tracemalloc
    ks = KondoSpectrum(_chain(n=8), 0, 1000.)
    assert len(cd._occupied_states(ks)) == ks.dim == 256
    eV = np.linspace(-5e-3, 5e-3, 61)
    fb = cd.FBuilder(1000., omega0=20e-3, Gamma0=5e-6)
    ref = _old_potential(ks, eV, -0.05, 0.2, Fb=fb)
    F = _CountF(fb)
    del F.odd_part
    monkeypatch.setattr(cd, "CHUNK_ENTRIES", 2e5)
    tracemalloc.start()
    got = cd.third_order_potential_dIdV(ks, eV, -0.05, 0.2, Fb=F)
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    assert peak < 20e6
    assert _rel(got, ref) < 1e-13


# item 4: the coefficients are built once per call


def test_levi_civita_coefficients_are_fully_symmetric():
    """What the hub blocks rest on: Im coeff[i,f,m] is the same under every
    permutation of (i,f,m), for a real and for a complex H."""
    for complex_h in (False, True):
        ks = KondoSpectrum(_chain(n=4, complex_h=complex_h), 0, 2.)
        R = _old_coefficients(ks, np.arange(ks.dim))
        for perm in ((1, 0, 2), (0, 2, 1), (2, 1, 0), (1, 2, 0), (2, 0, 1)):
            assert _rel(R.transpose(perm), R) < 1e-14


def _count_builds(monkeypatch):
    built = []
    levi = cd._levi_civita
    def counting(S, real, rows, cols, **kw):
        out = levi(S, real, rows, cols, **kw)
        built.append(out.size)
        return out
    monkeypatch.setattr(cd, "_levi_civita", counting)
    return built


@pytest.mark.parametrize("T,cplx", [(0.5, False), (2.0, False), (2.0, True),
                                    (1000., False)])
def test_hub_blocks_equal_the_single_pass(monkeypatch, T, cplx):
    """Blocks of one hub each (COEFF_BYTES=1): the numbers of the old single
    pass, and every coefficient built once, n_occ n_f dim for the occupied
    hubs and (n_f - n_occ) n_occ dim for the other final states."""
    ks = KondoSpectrum(_chain(complex_h=cplx), 0, T)
    occ = cd._occupied_states(ks)
    nocc, nf, dim = len(occ), len(cd._final_states(ks, occ, EV)), ks.dim
    assert nocc > 1
    fb = cd.FBuilder(T, omega0=20e-3, Gamma0=5e-6)
    ref = _old_third(ks, EV, -0.05, Fb=fb)
    monkeypatch.setattr(cd, "COEFF_BYTES", 1)
    built = _count_builds(monkeypatch)
    got = cd.third_order_kondo_dIdV(ks, EV, -0.05, Fb=fb)
    assert _rel(got, ref) < 1e-12
    assert sum(built) == nocc*nf*dim + (nf - nocc)*nocc*dim
    # the exact form (no odd_part) gives the same numbers
    F = _CountF(fb)
    del F.odd_part
    assert _rel(cd.third_order_kondo_dIdV(ks, EV, -0.05, Fb=F), ref) < 1e-12


def test_hub_blocks_build_once_whatever_the_number_of_biases(monkeypatch):
    """The saving: the old code built the blocks again for every group of
    biases its exchange accumulators held; the count no longer depends on
    the number of biases."""
    ks = KondoSpectrum(_chain(), 0, 1000.)
    dim = ks.dim
    monkeypatch.setattr(cd, "COEFF_BYTES", 16*4*dim*dim) # blocks of 4 hubs
    built = _count_builds(monkeypatch)
    fb = cd.FBuilder(1000., omega0=20e-3, Gamma0=5e-6)
    cd.third_order_kondo_dIdV(ks, EV, -0.05, Fb=fb)
    few = sum(built)
    built.clear()
    many = np.linspace(-5e-3, 5e-3, 301)
    cd.third_order_kondo_dIdV(ks, many, -0.05, Fb=fb)
    assert sum(built) == few == dim**3


# item 6: the series and the tails of Theta


def _theta_exact(x):
    # mpmath is not a dependency of dmrgpy, so only these tests skip
    mp = pytest.importorskip("mpmath")
    mp.mp.dps = 60
    x = mp.mpf(float(x))
    if abs(x) < 1:
        return mp.mpf(1)/2 + (mp.sinh(x) - x)/(2*(mp.cosh(x) - 1))
    # the same closed form without the cancellation of 1/2 against 1/2:
    # Theta(y) = e^y (e^y - 1 - y)/(e^y - 1)^2 for y < 0
    y = -abs(x)
    low = mp.exp(y)*(mp.exp(y) - 1 - y)/(mp.exp(y) - 1)**2
    return low if x < 0 else 1 - low


@pytest.mark.parametrize("x", [1e-7, 5e-5, 9.9e-5, -9.9e-5, -3e-5])
def test_theta_series_below_1e_4(x):
    """The cubic coefficient is -1/180: -1/360 was 2.7e-15 off at 9.9e-5."""
    got = float(cd.Theta(np.array([x]))[0])
    assert abs(got - float(_theta_exact(x))) < 2.5e-16


@pytest.mark.parametrize("x", [-40.5, -41., -60., -200.])
def test_theta_tail_above_40(x):
    """Theta(-|x|) -> (|x|-1) e^-|x|, to every digit (the +1 was 5% off)."""
    got = float(cd.Theta(np.array([x]))[0])
    ref = _theta_exact(x)
    assert abs(got - float(ref)) < 1e-14*float(ref)
    up = float(cd.Theta(np.array([-x]))[0])
    assert abs(up - float(_theta_exact(-x))) < 2.5e-16


@pytest.mark.parametrize("x,tol", [(2e-4, 1e-12), (1e-3, 2e-13), (0.1, 1e-15),
                                   (3., 1e-15), (39., 1e-15)])
def test_theta_closed_form_between(x, tol):
    """Left as it is: the closed form loses about 1e-16/|x| to cancellation
    just above the series."""
    for s in (1., -1.):
        got = float(cd.Theta(np.array([s*x]))[0])
        assert abs(got - float(_theta_exact(s*x))) < tol
