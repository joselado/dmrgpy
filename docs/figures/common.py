"""Helpers shared by every docs/figures/<id>/compute.py.

    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import common

    params = common.start(dict(model="...", n=16, maxm=[...]))  # prints them
    ...                                                           # compute
    common.save(__file__, params, energy=e, maxm=m)               # data.npz

save() writes data.npz next to compute.py, with the parameters as one JSON
string under the key "params", completed with the git commit, the host,
the dmrgpy path and the wall time. Read it back with common.load(path).
"""
import json
import os
import socket
import subprocess
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_T0 = time.time()


def git_commit():
    """The commit of this checkout, or the one d22.sh push recorded.

    The remote copy on D22 is not a git repository, so `d22.sh push`
    writes the local `git rev-parse HEAD` (and a "-dirty" mark) into
    docs/figures/COMMIT.
    """
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=HERE,
                             capture_output=True, text=True, check=True)
        commit = out.stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--", "src"],
                               cwd=HERE, capture_output=True, text=True).stdout
        return commit + ("-dirty" if dirty.strip() else "")
    except Exception:
        try:
            return open(os.path.join(HERE, "COMMIT")).read().strip()
        except OSError:
            return "unknown"


def start(params):
    """Print the parameters before the calculation starts; return them."""
    global _T0
    _T0 = time.time()
    print("parameters:")
    for k, v in params.items():
        print("  %s = %r" % (k, v))
    print(flush=True)
    return params


def save(script, params, **arrays):
    """Write data.npz next to `script`, with params stored as JSON."""
    import dmrgpy
    p = dict(params)
    p["git_commit"] = git_commit()
    p["host"] = socket.gethostname()
    p["dmrgpy_path"] = os.path.dirname(dmrgpy.__file__)
    p["wall_time_s"] = round(time.time() - _T0, 1)
    path = os.path.join(os.path.dirname(os.path.abspath(script)), "data.npz")
    np.savez(path, params=json.dumps(p), **arrays)
    print("wrote", path, "wall time %.1f s" % p["wall_time_s"])
    return path


def load(path):
    """Return (params dict, npz arrays) of a data.npz."""
    d = np.load(path)
    return json.loads(str(d["params"])), d
