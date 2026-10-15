"""Provenance of the user-guide figures.

A docs/figures/<id>/compute.py is plain dmrgpy and numpy: it ends with
np.savez("data.npz", ...) and knows nothing about this module. run.sh
calls stamp() after it, which adds the git commit, the host, the dmrgpy
path and the wall time to data.npz, as one JSON string under the key
"provenance". provenance(path) reads it back.
"""
import json
import os
import socket
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


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


def stamp(path, wall_time_s):
    """Add the provenance to the data.npz at `path`."""
    import dmrgpy
    with np.load(path) as f:
        arrays = {k: f[k] for k in f.files if k != "provenance"}
    p = dict(git_commit=git_commit(), host=socket.gethostname(),
             dmrgpy_path=os.path.dirname(dmrgpy.__file__),
             wall_time_s=float(wall_time_s))
    np.savez(path, provenance=json.dumps(p), **arrays)
    print("stamped", path, p)


def provenance(path):
    """The provenance dict of a data.npz."""
    with np.load(path) as f:
        return json.loads(str(f["provenance"]))


if __name__ == "__main__":
    stamp(sys.argv[1], sys.argv[2])
