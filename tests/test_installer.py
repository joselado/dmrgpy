"""The logic of install.py that needs no build: the compiler flags, the
version parser, the order of the candidates, the files it writes and its
command line. A real build is tested by install.py's own self-test."""
import os
import shutil
import sys
import types

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from installtk import blaslapack, buildflags, cppversion, install2  # noqa: E402
from installtk import requirements  # noqa: E402


# first lines of `--version`, as the compilers print them (the Apple and
# Red Hat ones are written from memory of their format)
@pytest.mark.parametrize("line, version", [
    ("g++ (Ubuntu 13.2.0-23ubuntu4) 13.2.0", (13, 2, 0)),
    ("g++-6 (Ubuntu 6.5.0-2ubuntu1~18.04) 6.5.0 20181026", (6, 5, 0)),
    ("g++ (GCC) 11.4.1 20231218 (Red Hat 11.4.1-3)", (11, 4, 1)),
    ("x86_64-conda-linux-gnu-g++ (conda-forge gcc 13.3.0-1) 13.3.0",
     (13, 3, 0)),
    ("Apple clang version 15.0.0 (clang-1500.3.9.4)", (15, 0, 0)),
    ("Ubuntu clang version 18.1.3 (1ubuntu1)", (18, 1, 3)),
    ("no numbers here", None),
])
def test_version_parser_reads_gcc_and_clang_lines(line, version):
    assert cppversion.parse_version(line) == version


def test_flags_on_x86_gcc_are_the_ones_the_builds_always_had():
    """Numbers must not move: on x86-64 with GCC the compiler command is
    the string install.py wrote into options.mk before the flags had a
    module of their own."""
    assert buildflags.cccom("g++", 3, "gcc", "x86_64") == \
        "g++ -m64 -std=c++17 -fconcepts -fPIC"
    assert buildflags.cccom("g++", 2, "gcc", "x86_64") == \
        "g++ -m64 -std=c++14 -fPIC"


def test_clang_gets_no_fconcepts_and_arm_gets_no_m64():
    assert buildflags.cccom("clang++", 3, "clang", "arm64") == \
        "clang++ -std=c++17 -Wno-gcc-compat -fPIC"
    assert "-m64" not in buildflags.cccom("g++", 3, "gcc", "aarch64")
    assert "-fconcepts" in buildflags.cccom("g++", 3, "gcc", "aarch64")
    with pytest.raises(ValueError):
        buildflags.std_flags(4)


def test_default_jobs_is_small_and_positive():
    assert 1 <= buildflags.default_jobs() <= 4
    assert buildflags.default_jobs(limit=1) == 1


@pytest.mark.skipif(shutil.which("g++") is None, reason="no g++ here")
@pytest.mark.parametrize("version", [2, 3])
def test_trial_compile_uses_the_standard_of_the_backend(version):
    family = cppversion.compiler_family("g++")
    flags = buildflags.std_flags(version, family).split()
    ok, out = cppversion.trial_compile("g++", std_flags=flags,
            source=buildflags.TRIAL_SOURCE[version])
    assert ok, out


@pytest.mark.skipif(shutil.which("g++") is None, reason="no g++ here")
def test_the_v3_trial_source_fails_under_the_old_cxx14_check():
    """What the check missed before: it compiled with -std=c++14 whatever
    backend was asked for."""
    ok, _ = cppversion.trial_compile("g++", std_flags=["-std=c++14"],
            source=buildflags.TRIAL_SOURCE[3])
    assert not ok


def test_explicit_gpp_is_the_only_compiler_candidate():
    got = list(requirements._compiler_candidates("my-g++", True, "python"))
    assert got == [("explicit --gpp", "my-g++")]


def test_conda_compiler_comes_before_the_system_one(monkeypatch):
    monkeypatch.setattr(cppversion, "find_conda_compiler",
            lambda python_exe=None: "/conda/bin/x-conda-g++")
    monkeypatch.setattr(requirements.shutil, "which",
            lambda name: "/usr/bin/"+name)
    got = list(requirements._compiler_candidates(None, True, "python"))
    assert [g[1] for g in got] == ["/conda/bin/x-conda-g++", "/usr/bin/g++"]
    got = list(requirements._compiler_candidates(None, False, "python"))
    assert [g[1] for g in got] == ["/usr/bin/g++"]


def test_blas_candidates_keep_their_order(monkeypatch):
    monkeypatch.setenv("CONDA_PREFIX", "/conda")
    monkeypatch.setattr(blaslapack, "_system_libdir_for",
            lambda libname, sys_gpp="g++": None)
    got = list(blaslapack.candidates("g++"))
    assert got[0]["libflags"] == "-lpthread -L/conda/lib -lblas -llapack"
    assert got[1]["libflags"] == "-lpthread -lblas -llapack"
    assert all(c["platform"] == "openblas" for c in got[2:])
    only = list(blaslapack.candidates("g++", openblas=True,
            openblas_libdir="/x/lib"))
    assert len(only) == 1 and only[0]["platform"] == "openblas"


_TEMPLATE = ("# a comment\nCCCOM=g++ -m64 -std=c++17 -fconcepts -fPIC\n"
             "PLATFORM=macos\nBLAS_LAPACK_LIBFLAGS=-framework Accelerate\n"
             "OPTIMIZATIONS=-O2 -DNDEBUG\n")


def test_options_mk_substitutes_three_lines_and_nothing_else():
    config = types.SimpleNamespace(platform="lapack",
            libflags="-lpthread -L/opt/a\\b -lblas -llapack", includeflags="")
    out = install2.options_mk(_TEMPLATE, "/c/g++ -m64 -std=c++14 -fPIC",
            config)
    assert out.split("\n") == ["# a comment",
            "CCCOM=/c/g++ -m64 -std=c++14 -fPIC", "PLATFORM=lapack",
            "BLAS_LAPACK_LIBFLAGS=-lpthread -L/opt/a\\b -lblas -llapack",
            "OPTIMIZATIONS=-O2 -DNDEBUG", ""]
    config.includeflags = "-I/x -DHAVE_LAPACK_CONFIG_H"
    out = install2.options_mk(_TEMPLATE, "g++", config)
    assert "\nBLAS_LAPACK_INCLUDEFLAGS=-I/x -DHAVE_LAPACK_CONFIG_H\n" in out


def test_build_config_mk_names_the_interpreter_and_the_compiler():
    text = install2.build_config_mk("g++ -m64 -std=c++14 -fPIC", "/p/python")
    assert "\nDMRGPY_PYTHON=/p/python\n" in text
    assert "\nDMRGPY_PYBIND_CCCOM=g++ -m64 -std=c++14 -fPIC\n" in text


def test_both_makefiles_read_build_config_and_stamp_the_record():
    for version in (2, 3):
        text = open(os.path.join(install2.backend_dir(version, ROOT),
                "Makefile")).read()
        assert "-include build_config.mk" in text
        assert "PYTHON=$(DMRGPY_PYTHON)" in text
        assert "PYBIND_CCCOM=$(DMRGPY_PYBIND_CCCOM)" in text
        assert "../buildinfo.py --stamp ." in text
        assert "$(wildcard *.h extra/*.h TDVP/*.h)" in text


def test_itensor_key_follows_the_configuration(monkeypatch):
    config = types.SimpleNamespace(platform="lapack", libflags="-lblas",
            includeflags="")
    monkeypatch.setattr(install2, "_git", lambda args, cwd:
            "abc" if args[0] == "rev-parse" else "")
    key = install2.itensor_key("ITensor", "g++ -O2", config)
    assert key and key == install2.itensor_key("ITensor", "g++ -O2", config)
    assert key != install2.itensor_key("ITensor", "clang++", config)
    config.libflags = "-lopenblas"
    assert key != install2.itensor_key("ITensor", "g++ -O2", config)
    # a compiler upgraded in place keeps its path, not its version
    config.libflags = "-lblas"
    assert key == install2.itensor_key("ITensor", "g++ -O2", config)
    config.compiler_version = (14, 1, 0)
    assert key != install2.itensor_key("ITensor", "g++ -O2", config)
    config.compiler_version = None
    # local changes in the vendored folder, or no git: never skip the build
    monkeypatch.setattr(install2, "_git", lambda args, cwd:
            "abc" if args[0] == "rev-parse" else " M itensor/x.cc")
    assert install2.itensor_key("ITensor", "g++ -O2", config) is None
    monkeypatch.setattr(install2, "_git", lambda args, cwd: None)
    assert install2.itensor_key("ITensor", "g++ -O2", config) is None


def test_clang_is_not_asked_for_cc1plus(monkeypatch):
    """Only GCC has a cc1plus; asking Clang for it rejected every Clang."""
    monkeypatch.setattr(requirements.shutil, "which", lambda name: name)
    monkeypatch.setattr(cppversion, "compiler_version", lambda gpp: (15, 0))
    monkeypatch.setattr(cppversion, "compiler_family", lambda gpp: "clang")
    monkeypatch.setattr(cppversion, "backend_missing", lambda gpp: True)
    seen = []
    monkeypatch.setattr(cppversion, "trial_compile",
            lambda gpp, std_flags=(), source=None, **kw:
            (seen.append(list(std_flags)) or (True, "")))
    assert requirements._try_compiler("clang++", True, (3,)) == (True, None)
    assert "-fconcepts" not in seen[0] and "-std=c++17" in seen[0]
    monkeypatch.setattr(cppversion, "compiler_family", lambda gpp: "gcc")
    ok, reason = requirements._try_compiler("g++", True, (3,))
    assert not ok and "cc1plus" in reason


def test_command_line_keeps_the_old_spellings():
    import install
    parser = install.build_parser()
    args = parser.parse_args(["--itensor_version=both", "--openblas",
            "--openblas_libdir=/a", "--openblas-includedir=/b",
            "--no-check-gpp", "--check-only"])
    assert args.itensor_version == "both" and args.doctor
    assert (args.openblas_libdir, args.openblas_includedir) == ("/a", "/b")
    assert args.check_gpp is False
    args = parser.parse_args([])
    assert args.link and args.selftest and not args.force
    # the shell rc file is edited only on request
    assert args.shellrc is False
    assert parser.parse_args(["--shellrc"]).shellrc is True
    assert parser.parse_args(["--no-shellrc"]).shellrc is False
    assert args.itensor_version == "3" and args.jobs is None
