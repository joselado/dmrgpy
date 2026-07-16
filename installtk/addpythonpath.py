import os
import sys
import sysconfig


def _site_packages():
    purelib = sysconfig.get_paths().get("purelib")
    if purelib and os.path.isdir(purelib):
        return purelib
    for p in sys.path:
        if "/site-packages" in p:
            return p
    return None


def addpath():
    mpath = os.path.dirname(os.path.realpath(__file__))
    sp = _site_packages()
    if sp is None:
        print("Could not find a site-packages directory to link into for "
              +sys.executable+".")
        print("Add this to your PYTHONPATH manually instead:")
        print("  "+mpath+"/../src")
        return

    target = os.path.join(sp, "dmrgpy")
    source = os.path.realpath(os.path.join(mpath, "..", "src", "dmrgpy"))

    if os.path.islink(target) or not os.path.exists(target):
        try:
            if os.path.islink(target) or os.path.exists(target):
                os.remove(target)
            os.symlink(source, target)
        except PermissionError:
            print("No permission to write to "+sp+".")
            print("Add this to your PYTHONPATH manually instead:")
            print("  "+mpath+"/../src")
            return
    else:
        print("Warning: "+target+" already exists and is not a symlink "
              "created by this installer -- leaving it untouched.")
        print("Add this to your PYTHONPATH manually instead:")
        print("  "+mpath+"/../src")
        return

    print("Linked dmrgpy into "+sp)
