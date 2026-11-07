"""z-averaged STM spectrum of a Kondo impurity on a superconductor, from
the discrete poles of the tangent-space Krylov (TASK) correlator.

One logarithmic Wilson mesh gives a few discrete levels per decade, so
above the gap the coherence peaks and the continuum of dI/dV are a comb of
poles. The cure of the numerical renormalization group (Bulla, Costi and
Pruschke, RMP 80, 395) is used here: several twisted meshes z (the twist
keyword of wilsonchain.star_levels), one chain and one DMRG + TASK solve
per z, and every pole broadened with the log-Gaussian kernel of
Weichselbaum and von Delft, PRL 99, 076402,

    K(x, x') = exp(-(ln(x/x')/b - b/4)^2) / (b x sqrt(pi)),   x, x' > 0,

which has unit area for any b and keeps a flat spectrum sampled on a
logarithmic mesh flat (the shift b/4 does that). Below omega0 it crosses
over to a Gaussian of width omega0, with the crossover weight
h = exp(-(ln(x'/omega0)/b)^2) for x' < omega0 read at the pole x' (the
"peak frequency" form of NRG Ljubljana), so every pole keeps exactly its
weight. The kernel is written from the paper; NRG Ljubljana's own code is
GPL-3 and is not copied.

The gap edge. A continuum pole at |w| > Delta is broadened in the
normal-state variable xi = sqrt(w^2 - Delta^2) and the density is mapped
back, A(w) = A_xi(xi) |w|/xi, the change of variables of the BCS density
of states. On a Wilson chain with uniform onsite pairing every Bogoliubov
level is exactly sqrt(e_k^2 + Delta^2) with e_k a level of the normal
chain, so this puts the log-Gaussian on the logarithmic mesh it was made
for, keeps every bit of weight outside the gap, and gives back the exact
|w|/sqrt(w^2-Delta^2) singularity from a smooth xi density.
The Gaussian part is reflected at xi = 0 so that it stays outside the gap.
With Delta = 0 the map is the identity and the kernel is the usual NRG
one on each sign of w. Poles with |w| < Delta are the subgap (YSR) poles:
they never enter the log kernel, are returned as a list, and are added to
the curve as Lorentzians of half-width delta_sub.

Normalization. The z-averaged first-site density of states of the
mean-energy mesh of wilsonchain.star_levels is rho/A_Lambda between
D Lambda^-1 and the lowest scale of the chain (A_Lambda =
ysr.lambda_factor), the same factor by which a coupling on the chain acts
weaker (ysr.lambda_factor's docstring). The spectrum of an operator
linear in f_0 is rho F(rho J, ...), so with J*A_Lambda on the chain the
continuum spectrum is A_Lambda times the chain one; zaverage_stm_spectrum
multiplies by it (continuum_norm=True). The top interval of a twisted
mean-energy mesh, [D Lambda^-z, D], is not corrected (the band-edge
artefact that the discretization of Zitko and Pruschke, PRB 79, 085106,
arXiv:0810.4765, removes; it is mesh="zp" of wilsonchain), so the
spectrum within a factor Lambda of D is not the continuum one: the scaled
curve carries A_Lambda/2 = 0.578 (Lambda=4) on each side against 0.5, the
bulk is right, and the excess sits within a factor Lambda of D and above
D (measured on the BdG chain, 20 sites, nz=1 and 4). On mesh="zp" the
z-averaged density is rho itself, the factor is 1 (kondo_chain's
fc.dos_factor, read here), and each side carries 0.500.

The gap edge on a short chain. The default omega0 of zaverage_stm_spectrum
is the last hopping of the chain, capped at the smallest xi of the
continuum poles when Delta != 0. So no continuum pole takes the Gaussian
part, the xi density vanishes at xi = 0 and the curve goes to zero at
|w| = Delta. With a larger omega0 given by hand the xi density is finite
at xi = 0 and the map imposes a 1/sqrt(|w| - Delta) divergence there:
within omega0 of the edge such a curve is not to be trusted, and a grid
point within roundoff of Delta picks up a huge value. A chain that does
not reach below Delta (fc.reaches_gap False) has no clean subgap window:
its lowest levels are strongly coupled to the impurity, and the YSR pole
can sit above Delta (measured: 1.12 Delta at 5 sites, J=0.6,
Delta=0.05), where it is broadened as continuum. zaverage_stm_spectrum
warns on such a chain. The sign of Delta is a gauge; only |Delta| is used.

Accuracy. The log-Gaussian is exact for a flat density and biased by
exp(p(p+1) b^2/4) where the density in xi goes locally as xi^p. At J = 0
the xi density is flat, and the z average of the BdG chain (Lambda=4, 20
sites) is the BCS density of states to 4e-6 at nz=4, b=0.4, 1.2e-3 at
nz=2, b=0.6 and 1.2e-2 at nz=1, b=1.0, from the resolution of the chain
up to D/Lambda. With a bound state inside the gap the density at the
impurity vanishes as sqrt(w - Delta) at the edge, p = 2, and there the
broadened curve is high by exp(3 b^2/2): 1.27 at b=0.4, 1.035 at b=0.15
(measured on a classical spin against the Shiba formula, 1.270, 1.034).
So b=0.4 with nz=4 is the choice for the continuum, and b of 0.15 with
nz of 16 when the shape within about 5% of Delta above the edge matters.

The TASK side. A Krylov space of nkry vectors resolves the lowest
excitations only to about their bandwidth over nkry^2, while the Wilson
levels near the gap edge, and near w = 0 in the normal state, are spaced
logarithmically far below that; the poles there are merged into a few
Ritz poles that carry the right total weight at shifted energies. At J = 0
on 20 sites (Delta = 1e-3, maxm=128) the density for xi in [Delta/10,
Delta] was off the BdG one of the same chain by up to 80% at nkry=100 and
by 4e-4 at nkry=300 (maxm=256 changed nothing). With rho J = 0.15 at
Delta = 1e-3, one mesh gave the same curve at nkry=300 and 600 from
w = 1.3 Delta up (to 3e-3), and a different one below 1.2 Delta. In the
normal state (Delta = 0, 21 sites, nkry=300) the lowest pole of each mesh
sat at 4e-5 to 6e-5, T_K/10, with the weight of everything below, so the
Kondo resonance below about T_K/3 is not resolved by this route.
nkry=300 is the default here; the cost grows linearly with it.

The ground state. Every tunneling vertex changes the fermion parity, and
near the singlet-doublet transition a plain DMRG solve can stop at a
mixture of the two parities, which gives an elastic weight and a split,
shifted subgap pole. parity=True (stm_poles, zaverage_stm_spectrum) takes
the poles on the ground state solved inside its parity sector instead
(ysrtk/paritystate.py, with the measurements), and parity="auto", the
default, does so only when two vevs on the plain state show a mixture. At J = 0 the
composite vertex leaves the one-site tangent space of the product ground
state, so TASK finds none of its spin-flip weight; stm_poles warns when
more than half of the weight is missing.

What was checked (tests/test_ysr_zaverage.py): the kernel's weight and
flatness, the gap map, the BdG poles of the chain against the TASK poles
of the many-body chain at J = 0, the z average of the BdG chain against
the BCS density of states and of a classical spin against the Shiba
formula, and the subgap pole against ysr.subgap_energies.
"""

import numpy as np

from .. import wilsonchain


def twists(nz):
    """The nz twists z_k = k/nz, k = 0..nz-1 (z = 0 is the plain mesh)."""
    if isinstance(nz, bool) or not isinstance(nz, (int, np.integer)) or nz < 1:
        raise ValueError("nz must be a positive integer, got %r" % (nz,))
    return np.arange(nz)/float(nz)


def log_gaussian(x, xp, b):
    """The log-Gaussian kernel K(x, x') of the module docstring, for x and
    x' positive (arrays broadcast)."""
    x = np.asarray(x, dtype=float)
    return (np.exp(-(np.log(x/xp)/b - b/4.0)**2)
            /(b*x*np.sqrt(np.pi)))


def _half_line_density(x, xp, wp, b, omega0):
    """Density on x > 0 of the poles xp >= 0 with weights wp."""
    x = np.asarray(x, dtype=float)
    out = np.zeros(len(x))
    if len(xp) == 0:
        return out
    xp = np.asarray(xp, dtype=float)
    wp = np.asarray(wp, dtype=float)
    with np.errstate(divide="ignore"):
        h = np.where(xp >= omega0, 1.0,
                     np.exp(-(np.log(np.maximum(xp, 1e-300)/omega0)/b)**2))
    good = x > 0.0
    xg = x[good][:, None]
    gauss = (np.exp(-((xg - xp)/omega0)**2) + np.exp(-((xg + xp)/omega0)**2))
    gauss = gauss/(omega0*np.sqrt(np.pi))
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        lg = np.where(xp > 0.0, log_gaussian(xg, np.maximum(xp, 1e-300), b), 0.0)
    K = h*lg + (1.0 - h)*gauss
    out[good] = K @ wp
    return out


def _default_omega0(omega, gap, t_last):
    """omega0 for one chain: its last hopping |t_last|, its energy
    resolution, capped at gap != 0 by the smallest xi of the continuum
    poles, so that no continuum pole takes the Gaussian part and the
    density vanishes at the gap edge (module docstring)."""
    w0 = abs(t_last)
    gap = abs(float(gap))
    if gap > 0.0:
        omega = np.abs(np.asarray(omega, dtype=float))
        xi = np.sqrt(omega[omega > gap]**2 - gap*gap)
        xi = xi[xi > 0.0]
        if len(xi):
            w0 = min(w0, float(np.min(xi)))
    return w0


def broaden_poles(es, omega, weight, b=0.4, omega0=None, gap=0.0,
                  delta_sub=None):
    """Density at the energies es of the poles omega with real weights
    weight, by the kernel of the module docstring.

    gap is the BCS gap Delta (only |gap| is used): poles with |omega| > gap are broadened in
    xi = sqrt(omega^2 - gap^2) and mapped back, poles with |omega| <= gap
    (subgap) are added as Lorentzians of half-width delta_sub (left out
    if delta_sub is None). omega0 is the energy, in the xi variable, below
    which the log-Gaussian becomes a Gaussian of width omega0 (default
    the smallest |xi| of the continuum poles). Returns (A, subgap), with
    subgap a dict of the subgap poles ("omega", "weight")."""
    es = np.asarray(es, dtype=float)
    omega = np.asarray(omega, dtype=float)
    weight = np.asarray(weight, dtype=float)
    if b <= 0.0:
        raise ValueError("b must be positive, got %r" % (b,))
    gap = abs(float(gap)) # the sign of Delta is a gauge
    sub = np.abs(omega) <= gap if gap > 0.0 else np.zeros(len(omega), bool)
    cont = ~sub
    xi = np.sqrt(np.maximum(omega[cont]**2 - gap*gap, 0.0))
    if omega0 is None:
        pos = xi[xi > 0.0]
        omega0 = np.min(pos) if len(pos) else 1e-12
    A = np.zeros(len(es))
    xe = np.sqrt(np.maximum(es*es - gap*gap, 0.0))
    oc, wc = omega[cont], weight[cont]
    zero = oc == 0.0 # only at gap = 0: half of such a pole goes to each side
    for sign in (1.0, -1.0):
        out = sign*es > gap
        if not np.any(out):
            continue
        on = sign*oc > 0.0
        xs = np.concatenate([xi[on], xi[zero]])
        ws = np.concatenate([wc[on], 0.5*wc[zero]])
        dens = _half_line_density(xe[out], xs, ws, b, omega0)
        if gap > 0.0: # back from xi to w, the BCS change of variables
            dens = dens*np.abs(es[out])/np.maximum(xe[out], 1e-300)
        A[out] = A[out] + dens
    subgap = dict(omega=omega[sub], weight=weight[sub])
    if delta_sub is not None and np.any(sub):
        lor = delta_sub/(np.pi*((es[:, None] - omega[sub][None, :])**2
                                + delta_sub**2))
        A = A + lor @ weight[sub]
    return A, subgap


def bdg_poles(eps, t, Delta, site=0):
    """Poles and weights of the spin-up local density of states of site
    site of a chain with onsite energies eps, hoppings t and onsite
    pairing Delta (a number, or one value per site), from its
    Bogoliubov-de Gennes matrix in the Nambu basis (f_up, f_dn^dag): one
    pole at every eigenvalue E, with weight |u_site(E)|^2, electron-like
    for E > 0 and hole-like for E < 0.
    The weights add up to 1."""
    n = len(eps)
    h = np.diag(eps) + np.diag(t, 1) + np.diag(t, -1)
    P = np.diag(np.broadcast_to(np.asarray(Delta, dtype=float), (n,)))
    H = np.block([[h, P], [P, -h]])
    E, V = np.linalg.eigh(H)
    return E, np.abs(V[site, :])**2


def bdg_chain_dos(es, n, Delta, Lambda=4.0, D=1.0, nz=4, b=0.4, omega0=None,
                  continuum_norm=True, mesh="mean"):
    """z-averaged, broadened local density of states per spin of the
    first site of the n-site Wilson chain of a BCS band, from the BdG
    poles: the single-particle reference of zaverage_stm_spectrum at J=0
    (where dI/dV is U^2/2 times this). With continuum_norm it is
    multiplied by A_Lambda on the default mesh (by 1 on mesh="zp"), so it
    stands for the continuum rho(w) = |w|/sqrt(w^2-Delta^2)/(2D)."""
    from ..ysr import lambda_factor
    A = np.zeros(len(es))
    for z in twists(nz):
        eps, t = wilsonchain.wilson_chain(n, Lambda=Lambda, band=(-D, D), z=z,
                                          mesh=mesh)
        w, wt = bdg_poles(eps, t, Delta)
        w0 = omega0
        if w0 is None:
            w0 = _default_omega0(w, Delta, t[-1] if len(t) else D)
        A = A + broaden_poles(es, w, wt, b=b, omega0=w0, gap=Delta)[0]
    A = A/nz
    if continuum_norm and mesh == "mean":
        A = A*lambda_factor(Lambda)
    return A


def stm_poles(fc, imp=None, site=None, U=0.0, nkry=300, operators=None,
              parity="auto", pad_bonds=None, expand_bonds=None,
              expand_rounds=1):
    """Discrete poles of dI/dV, A(w) of ysr.py's module docstring, from
    taskdmrg.get_poles: the electron part (X, X^dag) at +omega and the
    hole part (X^dag, X) at -omega, summed over the tip spin. imp and
    site are those of ysr.tunneling_operators (its defaults if None).
    operators replaces the tunneling vertices X_sigma^dag by a list of
    other creation-like operators (for instance [fc.Cdagup[1]] for the
    local density of states of spin up); U is then not used, and a
    nonzero U raises (fold it into the operators).

    Returns a dict with "omega", "weight" (real), "missing_weight", the
    sum of |missing_weight| of get_poles over the correlators (the weight
    the one-site tangent space dropped), and "elastic", the sum of the
    elastic weights |<X^dag><X>|. Every X changes the fermion parity, so
    the elastic weight vanishes on a ground state of definite parity and
    is left out of the poles; a nonzero one means that DMRG stopped at a
    mixture of the singlet and the doublet (measured: 5.5e-6 with
    <(S^z_tot)^2> = 1.8e-4 at maxm=48, nsweeps=15 on 10 sites), and a
    RuntimeWarning says so above 1e-6 of the total weight.

    With parity=True the poles are taken on the ground state solved
    inside its parity sector, paritystate.parity_ground_state(fc), so the
    mixture cannot appear; the result then also holds "parity", the info
    of parity_ground_state (sector energies, Sz_ground). fc itself is not
    solved then. With parity="auto", the default, fc is solved as with
    False and probed
    first (paritystate.auto_parity_state: two vevs per vertex, no
    correlator), and only when its elastic weight is above 1e-6 of the sum
    rule <X X^dag + X^dag X>, the threshold of the warning below, are the
    poles taken on the parity state as with True; the result then holds
    "parity_probe" (elastic, total, tripped, E0_plain), and "parity" only
    on a trip. parity takes False, True or "auto", nothing else. At
    J = 0 the ground state is a product state, the
    one-site tangent space misses X|GS> of the composite vertex, and TASK
    returns an empty spectrum with all of the weight missing; a
    RuntimeWarning says so when the missing weight is more than half of
    the total. pad_bonds (None, an int or "full") pads the bonds of the
    state before the TASK recursion, which enlarges the tangent space
    (taskdmrg.get_poles): "full" recovers the whole weight on a short
    chain at J = 0, and the spectrum of exact diagonalization on a
    promoted sector state of a degenerate pair, whose bonds are too small
    at the cutoff of the chain. On a long chain an int below the
    structural bound slowed the Krylov convergence of the subgap pole
    (pyitensor/task.py), so leave it None there. expand_bonds and
    expand_rounds go to taskdmrg.get_poles: they add to every bond at most
    expand_bonds directions taken from the part of X|GS> and H X|GS> the
    tangent space drops, which gives the spectrum of exact diagonalization
    on the promoted sector state, and with expand_rounds=2 on the product
    state at J = 0. They do not speed up the Krylov convergence of a pole
    on a long chain."""
    import warnings
    from .. import taskdmrg
    from ..ysr import tunneling_operators
    if operators is not None and U != 0.0:
        raise ValueError("stm_poles: U is used only to build the tunneling "
            "vertices; with operators given, fold U into them")
    from . import paritystate
    parity = paritystate.check_parity(parity)
    where = {k: v for k, v in (("imp", imp), ("site", site)) if v is not None}
    pinfo, probe = None, None
    if parity == "auto": # the plain solve, and the sectors only on a mixture
        ops = operators if operators is not None else \
            tunneling_operators(fc, U=U, **where)
        fc, pinfo, probe = paritystate.auto_parity_state(fc, ops)
    elif parity:
        fc, pinfo = paritystate.parity_ground_state(fc)
    if operators is None:
        operators = tunneling_operators(fc, U=U, **where)
    om, wt, missing, elastic = [], [], 0.0, 0.0
    for X in operators:
        Xd = X.get_dagger()
        for pair, sign in (((Xd, X), 1.0), ((X, Xd), -1.0)):
            p = taskdmrg.get_poles(fc, name=pair, nkry=nkry,
                                   pad_bonds=pad_bonds,
                                   expand_bonds=expand_bonds,
                                   expand_rounds=expand_rounds)
            om.append(sign*np.asarray(p["omega"], float))
            wt.append(np.real(p["weight"]))
            missing += abs(p["missing_weight"])
            elastic += abs(p["elastic"])
    om, wt = np.concatenate(om), np.concatenate(wt)
    if elastic > 1e-6*max(np.sum(np.abs(wt)), 1e-300):
        warnings.warn("stm_poles: elastic weight %.3g, but a tunneling "
            "vertex changes the fermion parity, so the ground state is a "
            "mixture of the two parities (DMRG not converged near the "
            "singlet-doublet degeneracy); pass parity=True to solve the "
            "ground state inside its parity sector, or raise nsweeps"
            % elastic,
            RuntimeWarning, stacklevel=2)
    if missing > 0.5*(missing + np.sum(np.abs(wt))):
        warnings.warn("stm_poles: the TASK poles miss %.3g of the weight "
            "(%.3g found). The part of X|GS> outside the one-site tangent "
            "space is dropped, which is all of it for the composite vertex "
            "on a product ground state (J = 0). Use expand_bonds=4 with "
            "expand_rounds=2, or pad_bonds=\"full\" or mode=\"ED\" on a "
            "short chain, or zaverage.bdg_chain_dos for the local density "
            "of states"
            % (missing, np.sum(np.abs(wt))), RuntimeWarning, stacklevel=2)
    out = dict(omega=om, weight=wt, missing_weight=missing, elastic=elastic)
    if pinfo is not None:
        out["parity"] = pinfo
    if probe is not None:
        out["parity_probe"] = probe
    return out


def zaverage_stm_spectrum(n, es, J=0.2, Delta=0.0, Lambda=4.0, D=1.0, nz=4,
                          b=0.4, omega0=None, U=0.0, delta_sub=None,
                          maxm=64, nsweeps=20, nkry=300, continuum_norm=True,
                          operators=None, parity_gap=False, parity="auto",
                          pad_bonds=None, expand_bonds=None,
                          expand_rounds=1, **kwargs):
    """dI/dV(eV) at the energies es, z-averaged over nz meshes
    (z = k/nz), each a ysr.kondo_chain(n, J=J, Delta=Delta, Lambda=Lambda,
    z=z) solved by DMRG with maxm and nsweeps, its poles from TASK
    (stm_poles) broadened by broaden_poles with width b, gap |Delta| and
    omega0 (default the last hopping of each chain, its energy
    resolution, capped at Delta != 0 by the smallest xi of its continuum
    poles, see the module docstring), and the nz curves averaged.

    J is the coupling on the chain, as in kondo_chain (J*A_Lambda stands
    for J in the continuum on the default mesh, J itself on mesh="zp").
    With continuum_norm the result is multiplied by fc.dos_factor
    (A_Lambda, or 1 on mesh="zp"), the normalization of the module
    docstring, so that at
    J = 0 it is U^2/2 times the continuum density of states per spin,
    rho(w) = |w|/sqrt(w^2-Delta^2)/(2D). Subgap poles are kept per z and
    added to the curve as Lorentzians of half-width delta_sub (left out
    if None). operators, a function of the chain that returns a list of
    creation-like operators, replaces the tunneling vertices (see
    stm_poles; U must then be 0). On a chain that does not reach below
    |Delta| (last hopping above |Delta|/10, kondo_chain's reaches_gap) a
    RuntimeWarning says that the subgap window is not clean (module
    docstring). With parity_gap, ysr.subgap_energies is solved on every
    chain as well, which gives the subgap energy to compare the subgap
    poles with. With parity, the poles of each mesh are taken on its
    ground state solved inside its parity sector
    (paritystate.parity_ground_state), so a parity mixture cannot appear;
    the plain solve is skipped, info["E0"] is the energy of that state,
    and E_ysr for parity_gap comes from the same sector solves. With
    parity="auto", the default, each mesh is solved plainly and probed
    (paritystate.auto_parity_state), and only a mesh whose plain state is
    a parity mixture takes the sector route of parity=True. pad_bonds,
    expand_bonds and expand_rounds go to stm_poles. Extra keywords go to kondo_chain.

    Returns (es, A, info): info["z"] the twists, and per z
    info["subgap"] (dict omega, weight of the subgap poles),
    info["missing_weight"], info["elastic"] (see stm_poles),
    info["omega0"], info["E0"], info["qp_gap"] (from kondo_chain) and
    info["reaches_gap"], info["A_z"] (the spectrum of that one mesh, with
    the same normalization as A), with parity_gap info["E_ysr"], and with
    parity info["parity"] (the sector energies, 2Sz -> energy). With
    parity="auto" the lists stay one entry per z: info["parity"] holds the
    sector energies of a mesh that tripped and None for one that did not,
    info["E0"] the sector E0 of a mesh that tripped and the plain E0
    otherwise, and info["parity_probe"] the probe of every mesh (elastic,
    total, tripped, E0_plain)."""
    import warnings
    from ..ysr import kondo_chain, subgap_energies, tunneling_operators
    from . import paritystate
    if operators is not None and U != 0.0:
        raise ValueError("zaverage_stm_spectrum: U is used only to build "
            "the tunneling vertices; with operators given, fold U into them")
    parity = paritystate.check_parity(parity)
    es = np.asarray(es, dtype=float)
    A = np.zeros(len(es))
    info = dict(z=twists(nz), subgap=[], missing_weight=[], elastic=[],
                omega0=[], E0=[], E_ysr=[], A_z=[], qp_gap=[],
                reaches_gap=[], parity=[])
    if parity == "auto":
        info["parity_probe"] = []
    for z in info["z"]:
        fc = kondo_chain(n, J=J, Delta=Delta, Lambda=Lambda, D=D, z=z, **kwargs)
        fc.maxm, fc.nsweeps = maxm, nsweeps
        pinfo = None
        if parity == "auto": # the plain solve, the sectors only on a mixture
            ops = None if operators is None else operators(fc)
            if ops is None:
                ops = tunneling_operators(fc, U=U)
            fc_gs, pinfo, probe = paritystate.auto_parity_state(fc, ops)
            info["parity_probe"].append(probe)
            info["E0"].append(probe["E0_plain"] if pinfo is None
                              else pinfo["E0"])
            info["parity"].append(None if pinfo is None
                                  else pinfo["sectors"])
        elif parity: # the state of definite parity, on its own copy of fc
            fc_gs, pinfo = paritystate.parity_ground_state(fc)
            info["E0"].append(pinfo["E0"])
            info["parity"].append(pinfo["sectors"])
        else:
            fc_gs = fc
            info["E0"].append(fc.gs_energy())
        ops = None if operators is None else operators(fc_gs)
        p = stm_poles(fc_gs, U=U, nkry=nkry, operators=ops,
                      pad_bonds=pad_bonds, expand_bonds=expand_bonds,
                      expand_rounds=expand_rounds)
        t_last = fc.t_last # the last hopping, on the mesh of the chain
        norm = fc.dos_factor if continuum_norm else 1.0
        w0 = omega0
        if w0 is None: # the energy resolution of this chain
            w0 = _default_omega0(p["omega"], Delta, t_last)
        Az, sub = broaden_poles(es, p["omega"], p["weight"], b=b, omega0=w0,
                                gap=Delta, delta_sub=delta_sub)
        A = A + Az
        info["A_z"].append(Az*norm)
        info["subgap"].append(sub)
        info["missing_weight"].append(p["missing_weight"])
        info["elastic"].append(p["elastic"])
        info["omega0"].append(w0)
        info["qp_gap"].append(getattr(fc, "qp_gap", None))
        # the criterion of kondo_chain's fc.reaches_gap
        info["reaches_gap"].append(bool(n > 1 and t_last < 0.1*abs(Delta)))
        if parity_gap:
            info["E_ysr"].append(pinfo["E_ysr"] if pinfo is not None
                                 else subgap_energies(fc)["E_ysr"])
    A = A*norm/nz # one mesh for all z, so one factor
    if Delta != 0.0 and not all(info["reaches_gap"]):
        gaps = [g for g in info["qp_gap"] if g is not None]
        warnings.warn("zaverage_stm_spectrum: the chain does not reach "
            "below Delta (last hopping above |Delta|/10, lowest bath "
            "quasiparticle up to %.4g |Delta|): its lowest levels are "
            "strongly coupled to the impurity, a YSR pole can sit above "
            "|Delta| and is then broadened as continuum, and the continuum "
            "near the edge is that of the chain; use a longer chain"
            % (max(gaps)/abs(Delta) if gaps else np.nan),
            RuntimeWarning, stacklevel=2)
    return es, A, info
