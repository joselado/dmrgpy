#!/usr/bin/env python3
"""Build, validate and publish the dmrgpy PyPI package.

Ships only the pure-Python part of dmrgpy (ED backend + the pyitensor
DMRG/TDVP backend) -- see the "Packaging / PyPI" section of CLAUDE.md for
what's deliberately excluded and why. Run from the repo root:

    python release_pypi.py                     # build + validate + upload to TestPyPI
    python release_pypi.py --repository pypi    # the real thing
    python release_pypi.py --no-upload          # build + validate only, upload nothing
    python release_pypi.py --skip-tests          # skip `pytest tests` (faster iteration)

Uploads always ask for interactive confirmation first, unless --yes is
given. A version can never be re-uploaded to (Test)PyPI once published --
bump pyproject.toml's [project].version for the next attempt.
"""
import argparse
import shutil
import subprocess
import sys
import tempfile
import venv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DIST = ROOT / "dist"
BUILD = ROOT / "build"


def run(cmd, **kwargs):
    print("+ " + " ".join(str(c) for c in cmd))
    subprocess.run(cmd, check=True, cwd=ROOT, **kwargs)


def read_version():
    import tomllib
    with open(ROOT / "pyproject.toml", "rb") as f:
        return tomllib.load(f)["project"]["version"]


def check_git_clean():
    out = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=no"],
        cwd=ROOT, check=True, capture_output=True, text=True,
    ).stdout
    if out.strip():
        sys.exit(
            "ERROR: tracked files have uncommitted changes -- commit or "
            "stash before releasing:\n" + out
        )


def check_tag_available(version):
    tag = f"v{version}"
    existing = subprocess.run(
        ["git", "tag", "-l", tag], cwd=ROOT, check=True, capture_output=True, text=True,
    ).stdout.strip()
    if existing:
        sys.exit(
            f"ERROR: tag {tag} already exists locally -- version {version} was "
            "already released. Bump [project].version in pyproject.toml first."
        )


def run_tests():
    run([sys.executable, "-m", "pytest", "tests"])


def clean_artifacts():
    for d in (DIST, BUILD):
        if d.exists():
            shutil.rmtree(d)
    for egg_info in ROOT.glob("src/*.egg-info"):
        shutil.rmtree(egg_info)


def build_package():
    run([sys.executable, "-m", "build"])


def twine_check():
    run([sys.executable, "-m", "twine", "check", *sorted(str(p) for p in DIST.glob("*"))])


def smoke_test_wheel():
    """Install the built wheel into a throwaway venv (not from src/, so
    missing package-data or import bugs actually surface) and cross-check
    the pure-Python DMRG backend against ED on a small Heisenberg chain."""
    wheels = sorted(DIST.glob("*.whl"))
    if not wheels:
        sys.exit("ERROR: no wheel found in dist/ to smoke-test")
    check_script = """
import dmrgpy.spinchain as spinchain

n = 4
spins = ["S=1/2"] * n

def heisenberg(sc):
    h = 0
    for i in range(n - 1):
        h = h - sc.Sz[i]*sc.Sz[i+1] - sc.Sx[i]*sc.Sx[i+1] - sc.Sy[i]*sc.Sy[i+1]
    sc.set_hamiltonian(h)
    return sc

c_dmrg = heisenberg(spinchain.Spin_Chain(spins))
c_dmrg.setup_python()
e_dmrg = c_dmrg.gs_energy(mode="DMRG")

e_ed = heisenberg(spinchain.Spin_Chain(spins)).gs_energy(mode="ED")

assert abs(e_dmrg - e_ed) < 1e-6, (e_dmrg, e_ed)
print(f"smoke test OK: pyitensor DMRG={e_dmrg!r} ED={e_ed!r}")
"""
    with tempfile.TemporaryDirectory() as tmp:
        venv_dir = Path(tmp) / "smoke-venv"
        venv.create(venv_dir, with_pip=True)
        pip = venv_dir / "bin" / "pip"
        python = venv_dir / "bin" / "python"
        run([str(pip), "install", "--quiet", str(wheels[0])])
        run([str(python), "-c", check_script])


def upload(repository):
    args = [sys.executable, "-m", "twine", "upload"]
    if repository != "pypi":
        args += ["--repository", repository]
    args += sorted(str(p) for p in DIST.glob("*"))
    run(args)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repository", choices=["testpypi", "pypi"], default="testpypi",
                         help="where to upload (default: testpypi)")
    parser.add_argument("--no-upload", action="store_true", help="build and validate only, skip the upload step")
    parser.add_argument("--skip-tests", action="store_true", help="skip running `pytest tests` first")
    parser.add_argument("--yes", action="store_true", help="skip the interactive upload confirmation")
    args = parser.parse_args()

    version = read_version()
    print(f"dmrgpy {version}")

    check_git_clean()
    check_tag_available(version)

    if not args.skip_tests:
        run_tests()

    clean_artifacts()
    build_package()
    twine_check()
    smoke_test_wheel()

    if args.no_upload:
        print(f"\nBuild OK: dist/dmrgpy-{version}* is ready. Not uploading (--no-upload).")
        return

    if not args.yes:
        reply = input(f"\nUpload dist/dmrgpy-{version}* to {args.repository}? [y/N] ")
        if reply.strip().lower() != "y":
            print("Aborted before upload.")
            return

    upload(args.repository)

    print(f"\nUploaded to {args.repository}.")
    if args.repository == "testpypi":
        print(
            "Verify with:\n"
            "  python -m venv /tmp/dmrgpy-testpypi-check && "
            "/tmp/dmrgpy-testpypi-check/bin/pip install "
            "--index-url https://test.pypi.org/simple/ "
            "--extra-index-url https://pypi.org/simple/ dmrgpy\n"
            "Once satisfied, publish for real with:\n"
            "  python release_pypi.py --repository pypi\n"
            f"Then tag the release:\n"
            f"  git tag v{version} && git push origin v{version}"
        )
    else:
        print(
            f"Tag the release:\n  git tag v{version} && git push origin v{version}\n"
            f"Remember: v{version} can never be re-uploaded to PyPI -- the next "
            "release needs a version bump."
        )


if __name__ == "__main__":
    main()
