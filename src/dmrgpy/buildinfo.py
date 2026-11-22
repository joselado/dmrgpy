"""Record of how a compiled C++ extension was built.

Each backend folder (mpscpp2/, mpscpp3/) holds a `build_info.json` next to
its `_dmrgcpp*.so`. `install.py` writes the configuration of the build
into it (compiler, flags, BLAS/LAPACK, interpreter, loaded modules,
commit), and the `pybind` target of the Makefile stamps it with a hash of
the sources the extension was compiled from. `cppext.py` reads it to say
why an extension does not load and whether it is older than its sources.

This file has no import from the rest of dmrgpy, so that the Makefile can
run it as a script:

    python ../buildinfo.py --stamp .
"""
import glob
import hashlib
import json
import os
import sys
import time

INFO_NAME = "build_info.json"

# the sources of the extension that this project owns; the vendored ITensor
# under ITensor/ is a separate static library and is not hashed
SOURCE_GLOBS = ("bindings.cc", "*.h", "extra/*.h", "TDVP/*.h")


def source_files(dirpath):
    """Sorted paths, relative to `dirpath`, of the extension's sources."""
    out = set()
    for pattern in SOURCE_GLOBS:
        for path in glob.glob(os.path.join(dirpath, pattern)):
            if os.path.isfile(path):
                out.add(os.path.relpath(path, dirpath))
    return sorted(out)


def sources_hash(dirpath):
    """sha256 of the names and contents of the extension's sources, or None
    when `dirpath` holds no sources (a pip install ships none)."""
    files = source_files(dirpath)
    if "bindings.cc" not in files:
        return None
    h = hashlib.sha256()
    for name in files:
        h.update(name.replace(os.sep, "/").encode())
        h.update(b"\0")
        with open(os.path.join(dirpath, name), "rb") as f:
            h.update(f.read())
        h.update(b"\0")
    return h.hexdigest()


def info_path(dirpath):
    return os.path.join(dirpath, INFO_NAME)


def read(dirpath):
    """The build record of `dirpath` as a dict, or None when there is none
    or it cannot be parsed."""
    try:
        with open(info_path(dirpath)) as f:
            info = json.load(f)
    except (OSError, ValueError):
        return None
    return info if isinstance(info, dict) else None


def write(dirpath, info):
    tmp = info_path(dirpath)+".tmp"
    with open(tmp, "w") as f:
        json.dump(info, f, indent=1, sort_keys=True)
        f.write("\n")
    os.replace(tmp, info_path(dirpath))


def stamp(dirpath, **fields):
    """Record that the extension in `dirpath` was just built from the
    sources that are there now. Keeps every other field of the record."""
    info = read(dirpath) or {}
    info.update(fields)
    info["sources_hash"] = sources_hash(dirpath)
    info["built"] = time.strftime("%Y-%m-%d %H:%M:%S")
    write(dirpath, info)
    return info


def is_stale(dirpath, info=None):
    """True when the record says the extension was built from other
    sources than the ones in `dirpath` now, False when they agree, and None
    when this cannot be decided (no record, no hash, or no sources)."""
    if info is None:
        info = read(dirpath)
    if not info or not info.get("sources_hash"):
        return None
    current = sources_hash(dirpath)
    if current is None:
        return None
    return current != info["sources_hash"]


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--stamp":
        stamp(os.path.abspath(sys.argv[2]))
    else:
        sys.exit("usage: python buildinfo.py --stamp <backend folder>")
