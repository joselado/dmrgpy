#!/usr/bin/python3
"""Set up the Julia backend of dmrgpy (itensor_version="julia_live").

It installs juliacall into the current Python interpreter if it is
missing, starts the Julia session once, which makes juliacall provision
Julia and the packages of src/dmrgpy/juliapkg.json (ITensors, ITensorMPS,
KrylovKit, VectorInterface), and links dmrgpy into the interpreter's
site-packages (--no-link skips that). The first run can take several
minutes: Julia downloads and precompiles its packages.
"""
import argparse
import importlib.util
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.realpath(__file__))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--link", default=True,
            action=argparse.BooleanOptionalAction,
            help="Link dmrgpy into the site-packages of this interpreter "
                 "(default: on).")
    args = parser.parse_args(argv)

    if importlib.util.find_spec("juliacall") is None:
        print("juliacall not found, installing it into the current "
              "environment ("+sys.executable+")...")
        proc = subprocess.run([sys.executable, "-m", "pip", "install",
                "juliacall"])
        if proc.returncode != 0:
            print("Could not install juliacall. Install it by hand with:\n"
                  "  "+sys.executable+" -m pip install juliacall")
            return 1

    print("Setting up the Julia backend (juliacall). This may take a few "
          "minutes the first time.")
    sys.path.insert(0, os.path.join(ROOT, "src"))
    ok = True
    try:
        import dmrgpy.mpsjulialive.juliasession # noqa: F401
        print("Julia backend ready")
    except Exception as e:
        ok = False
        print("Could not set up the Julia backend:", e)
        print("setup_julia() will not be usable until this is resolved")

    if args.link:
        from installtk import addpythonpath
        addpythonpath.addpath() # add the library to the Python path
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
