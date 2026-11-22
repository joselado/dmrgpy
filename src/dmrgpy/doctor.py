"""Say which DMRG backends can be used here, and why one cannot.

    python -m dmrgpy.doctor               # the state of every backend
    python -m dmrgpy.doctor --selftest    # also solve a small chain on each
                                          # compiled backend, against ED

`install.py` runs the self-test in a fresh process at the end of a build,
so that an extension that compiles and does not load, or loads and gives a
wrong energy, is found there and not later in a calculation.
"""
import argparse
import importlib.util
import os
import sys

from . import cppext

SELFTEST_SITES = 6
SELFTEST_TOL = 1e-6


def _build_lines(build):
    """The lines that describe a build record."""
    if not build:
        return ["no build record (built by hand, or before install.py "
                "wrote one)"]
    out = []
    first = "built "+str(build.get("built", "at an unknown time"))
    if build.get("commit"):
        first += ", commit "+str(build["commit"])[:10]
    out.append(first)
    if build.get("compiler"):
        out.append("compiler  "+str(build["compiler"])+"  ("
                +str(build.get("compiler_version", "?"))+")")
    if build.get("cccom"):
        out.append("flags     "+str(build["cccom"]))
    if build.get("blas"):
        out.append("BLAS      "+str(build["blas"])+"  "
                +str(build.get("libflags", "")))
    if build.get("python"):
        out.append("Python    "+str(build["python"]))
    if build.get("modules"):
        out.append("modules   "+", ".join(build["modules"]))
    if build.get("conda_prefix"):
        out.append("conda     "+str(build["conda_prefix"]))
    return out


def report():
    """The state of every backend, as text."""
    pad = " "*17
    lines = ["dmrgpy     "+os.path.dirname(os.path.abspath(__file__)),
             "Python     "+sys.executable+"  ("+sys.version.split()[0]+")",
             "threads    "+", ".join(name+"="+os.environ.get(name, "unset")
                     for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS")),
             ""]
    for version in (3, 2):
        st = cppext.status(version)
        name = ("ITensor v%s"%version).ljust(15)
        if st["available"]:
            lines.append(name+"  OK  "+str(st["file"]))
            if st["stale"]:
                lines.append(pad+"OLDER THAN ITS SOURCES: rebuild with "
                        "python install.py --extension-only")
            lines += [pad+l for l in _build_lines(st["build"])]
        else:
            lines.append(name+"  NOT AVAILABLE")
            lines.append(pad+cppext.explain(version))
            if st["state"] != cppext.NOT_COMPILED:
                lines += [pad+l for l in _build_lines(st["build"])]
    lines.append("python".ljust(15)+"  OK  (pure Python, always available)"
            if cppext.available("python")
            else "python".ljust(15)+"  NOT AVAILABLE  "+cppext.explain("python"))
    lines.append("julia_live".ljust(15)+("  juliacall is installed (the Julia "
            "session itself is not started by this check)"
            if importlib.util.find_spec("juliacall") is not None
            else "  NOT AVAILABLE  juliacall is not installed "
                 "(python install_julia.py)"))
    lines.append("")
    lines.append("A chain with no itensor_version runs on: "
            +repr(cppext.default_backend()))
    return "\n".join(lines)


def selftest(version):
    """Solve a small Heisenberg chain on the compiled backend `version` and
    on ED. Returns (ok, message)."""
    if not cppext.available(version):
        return False, cppext.explain(version)
    from . import mode, spinchain
    n = SELFTEST_SITES
    sc = spinchain.Spin_Chain(["S=1/2"]*n, itensor_version=version)
    sc.set_hamiltonian(sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
            + sc.Sz[i]*sc.Sz[i+1] for i in range(n-1)))
    if mode.get_mode(sc, mode="DMRG") != "DMRG":
        return False, "the chain did not run on DMRG"
    e_dmrg = sc.gs_energy(mode="DMRG")
    e_ed = sc.gs_energy(mode="ED")
    err = abs(e_dmrg-e_ed)
    msg = ("%d-site Heisenberg chain: E(DMRG) = %.10f, E(ED) = %.10f, "
           "difference %.1e"%(n, e_dmrg, e_ed, err))
    return err < SELFTEST_TOL, msg


def main(argv=None):
    parser = argparse.ArgumentParser(prog="python -m dmrgpy.doctor",
            description="Say which DMRG backends can be used here, and why "
                        "one cannot.")
    parser.add_argument("--selftest", nargs="*", type=int, choices=[2, 3],
            default=None, metavar="VERSION",
            help="solve a small chain on these compiled backends (2, 3) "
                 "and compare with exact diagonalization; with no version, "
                 "on every backend that is available")
    args = parser.parse_args(argv)
    print(report())
    if args.selftest is None:
        return 0
    versions = args.selftest or [v for v in (3, 2) if cppext.available(v)]
    if not versions:
        print("\nSelf-test: no compiled backend is available.")
        return 1
    failed = False
    print("")
    for version in versions:
        ok, msg = selftest(version)
        print("Self-test, ITensor v%s: %s  %s"%(version,
                "PASSED" if ok else "FAILED", msg))
        failed = failed or not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
