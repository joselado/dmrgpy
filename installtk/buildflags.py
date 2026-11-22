"""The compiler flags of each backend, in one place.

Both the requirement check (which trial-compiles with them) and the build
(which writes them into ITensor's options.mk and into the extension's
build_config.mk) read them from here, so that what is checked is what is
built.
"""
import os
import platform

VERSIONS = (2, 3)

# A small program that needs the language standard of each backend, so
# that a compiler too old for it fails in the check and not minutes into
# the ITensor build.
TRIAL_SOURCE = {
    2: "#include <memory>\n"
       "int main(){ auto p = std::make_unique<int>(0); return *p; }\n",
    3: "#include <optional>\n#include <utility>\n"
       "template<class T> int f(T x){ if constexpr (sizeof(T) > 1) "
       "return 0; else return x; }\n"
       "int main(){ std::optional<int> o = 0; "
       "auto [a, b] = std::pair<int,int>(*o, 0); return f(a + b); }\n",
}


def std_flags(version, family="gcc"):
    """Language-standard flags. ITensor v2 needs C++14. ITensor v3 needs
    C++17, with the concepts TS extension on GCC (-fconcepts); Clang has no
    such flag, and ITensor's own options.mk.sample builds there with
    -Wno-gcc-compat instead."""
    if version == 2:
        return "-std=c++14"
    if version == 3:
        if family == "clang":
            return "-std=c++17 -Wno-gcc-compat"
        return "-std=c++17 -fconcepts"
    raise ValueError("Unknown ITensor backend version: "+str(version))


def arch_flags(machine=None):
    """-m64 on x86-64, where the builds have always carried it, and nothing
    elsewhere: GCC on ARM does not know the flag."""
    machine = (machine or platform.machine()).lower()
    return "-m64" if machine in ("x86_64", "amd64") else ""


def cccom(gpp, version, family="gcc", machine=None):
    """The full compiler command of a backend: what ITensor's options.mk
    calls CCCOM, also used for the pybind11 extension. -fPIC is needed
    throughout, or the final link of the shared object fails with a TLS
    relocation error."""
    parts = [gpp, arch_flags(machine), std_flags(version, family), "-fPIC"]
    return " ".join(p for p in parts if p)


def default_jobs(limit=4):
    """Parallel compile jobs when none are asked for: the cores this
    process may use, at most `limit`. Small on purpose: an ITensor v3
    translation unit takes a lot of memory, and the build often runs on a
    shared login node."""
    try:
        n = len(os.sched_getaffinity(0))
    except AttributeError: # not on Linux
        n = os.cpu_count() or 1
    return max(1, min(limit, n))
