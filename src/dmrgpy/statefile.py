"""Wavefunctions on disk: one .npz file per state.

A state is written as plain NumPy arrays plus a JSON header, so that a
later process (a second cluster job, a notebook) reads it back onto a
chain built the same way:

    sc.get_gs().save("gs")          # or sc.save_state(wf,"gs")
    ...
    wf = sc.load_state("gs")        # a new process, the same chain
    sc.set_gs(wf)

An MPS is stored as its site tensors, (left link, site, right link) arrays
with length-1 links at the open ends, exactly as the session exports them
(mps_site_arrays on itensor_version=3 and "python", mpsjulialive/arrays.jl
on julia_live): no change of gauge and no change of norm. An ED state is
stored as its dense vector. The header holds the format version, the kind
of state, the site list of the chain, the conserved sector and the backend
that wrote the file.

What a file is tied to:

- the chain: the number of sites and the site list must be the ones of the
  chain that loads it, and for an ED state the conserved sector too, since
  its vector has the dimension of the sector;
- not the backend: the site arrays are in the basis of the site types,
  which itensor_version=3, "python" and julia_live share
  (tests/test_state_file.py loads every file on every backend);
- not the bond dimension: the state comes back with the bonds it has
  (the importers drop only exactly-zero singular values).

A state computed inside a conserved sector is written as its dense form
(the QN blocks scattered into dense arrays, what promote_mps() gives), and
an MPS loads on a chain with no sector set: the three array importers
build dense links only. itensor_version=2 has no array bridge and raises.
"""

import json
import os

import numpy as np

FORMAT = "dmrgpy-state"
VERSION = 1


def filename(name):
    """The file a state name stands for: `name` itself when it ends with
    .npz, name+".npz" otherwise (np.savez adds the suffix on its own, and
    np.load does not)."""
    name = os.fspath(name)
    return name if name.endswith(".npz") else name+".npz"


def _number(x):
    x = complex(x)
    if x.imag != 0.0:
        raise ValueError("a conserved charge must be real (got %r)" % (x,))
    x = x.real
    return int(x) if float(x).is_integer() else float(x)


def _sector(obj):
    """The conserved sector of a chain or of an ED object, as a plain
    dict, or None"""
    sec = getattr(obj, "conserved_sector", None)
    if not sec: return None
    return {str(k): _number(sec[k]) for k in sec}


def _site_labels(chain):
    sites = getattr(chain, "sites", None)
    if sites is None: return None
    return [str(s) for s in sites]


def _mps_arrays(chain, wf):
    """The site arrays of an MPS of `chain`, on the host, complex128. In
    a conserved sector too: every exporter makes the QN blocks dense."""
    from .mpsalgebra import _is_julia
    from .mpsalgebratk.hadamard import _host
    if _is_julia(wf):
        from .mpsjulialive.juliasession import Main
        return _host([np.asarray(a) for a in Main.mps_site_arrays(wf.jlmps)])
    session = getattr(chain, "_session", None)
    if (wf.cpp_handle is None or session is None
            or not hasattr(session, "mps_site_arrays")):
        raise NotImplementedError("saving an MPS needs itensor_version=3, "
                "'python' or 'julia_live' (got itensor_version=%r)"
                % (getattr(chain, "itensor_version", None),))
    return _host(session.mps_site_arrays(wf.cpp_handle))


def save_state(wf, name, chain=None):
    """Write the wavefunction `wf` to the file `name` (.npz is added when
    missing) and return the file name. `chain` is the chain of the state;
    an MPS knows its own, an ED state does not, and without it the file of
    an ED state has no site list to be checked against at load."""
    from .mpsalgebra import wavefunction_mode
    kind = wavefunction_mode(wf)
    meta = dict(format=FORMAT, version=VERSION)
    if kind == "ED":
        v = np.ascontiguousarray(np.asarray(wf.v), dtype=complex).reshape(-1)
        meta.update(kind="ED", dim=int(v.size),
                    sector=_sector(chain if chain is not None else wf.MBO),
                    sites=_site_labels(chain) if chain is not None else None,
                    backend="ED")
        arrays = {"v": v}
    else:
        if chain is None: chain = wf.MBO
        if chain is None:
            raise RuntimeError("the MPS has no chain attached; save it with "
                               "chain.save_state(wf,name)")
        arrs = _mps_arrays(chain, wf)
        meta.update(kind="MPS", ns=len(arrs), sites=_site_labels(chain),
                    sector=_sector(chain),
                    backend=str(chain.itensor_version))
        arrays = {"A%d" % i: a for i, a in enumerate(arrs)}
    meta["chain"] = type(chain).__name__ if chain is not None else None
    fname = filename(name)
    # through a temporary file in the same folder, so that a job killed
    # while it writes leaves the previous file and not half of a new one
    tmp = fname+".tmp%d" % os.getpid()
    try:
        with open(tmp, "wb") as f:
            np.savez(f, meta=np.array(json.dumps(meta)), **arrays)
        os.replace(tmp, fname)
    finally:
        if os.path.exists(tmp): os.remove(tmp)
    return fname


def read_header(name):
    """The header of a state file, as a dict"""
    with np.load(filename(name), allow_pickle=False) as f:
        return _header(f, filename(name))


def _header(f, fname):
    if "meta" not in f.files:
        raise ValueError("%s is not a dmrgpy state file" % fname)
    meta = json.loads(f["meta"].item())
    if meta.get("format") != FORMAT:
        raise ValueError("%s is not a dmrgpy state file" % fname)
    if int(meta.get("version", 0)) > VERSION:
        raise ValueError("%s has format version %s; this dmrgpy reads up "
                "to version %d" % (fname, meta.get("version"), VERSION))
    return meta


def _check_sites(chain, meta, fname):
    here = _site_labels(chain)
    there = meta.get("sites")
    if there is not None and here is not None and list(there) != here:
        raise ValueError("%s was written on a chain with sites %s; this "
                "chain has sites %s" % (fname, _short(there), _short(here)))


def _short(labels):
    labels = list(labels)
    if len(labels) <= 12: return "["+", ".join(labels)+"]"
    return "[%s, ... (%d sites)]" % (", ".join(labels[:6]), len(labels))


def _same_sector(a, b):
    a, b = a or {}, b or {}
    return (set(a) == set(b)
            and all(float(a[k]) == float(b[k]) for k in a))


def load_state(chain, name):
    """The wavefunction of the file `name`, as a state of `chain`: an MPS
    of the chain's backend, or an ED state, whichever the file holds. The
    state is returned as it was written, norm included, and nothing on the
    chain changes; set_gs(wf) or gs_energy(wf0=wf) make it the chain's."""
    fname = filename(name)
    with np.load(fname, allow_pickle=False) as f:
        meta = _header(f, fname)
        kind = meta.get("kind")
        if kind == "ED":
            data = np.array(f["v"], dtype=complex)
        elif kind == "MPS":
            data = [np.array(f["A%d" % i], dtype=complex)
                    for i in range(int(meta["ns"]))]
        else:
            raise ValueError("%s holds a state of unknown kind %r"
                             % (fname, kind))
    _check_sites(chain, meta, fname)
    if kind == "ED": return _load_ed(chain, data, meta, fname)
    return _load_mps(chain, data, meta, fname)


def _load_ed(chain, v, meta, fname):
    from .edtk.edchain import State
    if not _same_sector(meta.get("sector"), _sector(chain)):
        raise ValueError("%s holds an ED state of the sector %s; this "
                "chain has %s. An ED vector has the dimension of its "
                "sector, so set the same one with set_conserved_sector()"
                % (fname, meta.get("sector"), _sector(chain)))
    # the ED objects of the chain types have no common attribute for the
    # dimension, and every one of them can draw a random state
    try:
        ref = chain.get_ED_obj().random_state()
    except AttributeError:
        if getattr(chain, "hamiltonian", None) is not None: raise
        raise RuntimeError("%s holds an ED state, and the ED object of "
                "this chain is built from its Hamiltonian: call "
                "set_hamiltonian() before load_state()" % fname) from None
    if np.asarray(ref.v).size != v.size:
        raise ValueError("%s holds an ED vector of dimension %d; the ED "
                "space of this chain has dimension %d"
                % (fname, v.size, np.asarray(ref.v).size))
    return State(v, ref.MBO)


def _load_mps(chain, arrays, meta, fname):
    if len(arrays) != chain.ns:
        raise ValueError("%s holds an MPS of %d sites; this chain has %d"
                         % (fname, len(arrays), chain.ns))
    if getattr(chain, "conserved_sector", None):
        raise ValueError("conserved-sector mode (%s): an MPS loads on dense "
                "sites only. Load %s on the chain with the sector cleared "
                "(set_conserved_sector()); a state written inside a sector "
                "is stored as its dense form, the one promote_mps() gives"
                % (chain.conserved_sector, fname))
    from .gaussianmps import _bridge
    try:
        bridge = _bridge(chain)
    except NotImplementedError:
        raise NotImplementedError("loading an MPS needs itensor_version=3, "
                "'python' or 'julia_live' (got itensor_version=%r)"
                % (chain.itensor_version,)) from None
    return bridge.to_mps(arrays)
