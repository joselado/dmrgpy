"""Sweeps: per-sweep DMRG/TDVP parameters (bond dimension cap, truncation
cutoff, density-matrix noise, local-solver iteration count).

Real ITensor's Sweeps::maxdim()/cutoff()/noise() return a reference-like
proxy so C++ code can write `sweeps.maxdim() = 30` to broadcast one value
to every sweep, or `sweeps.maxdim(i)` to read one sweep's value. Python has
no such call-then-assign idiom, so this exposes maxdim/cutoff/noise/niter
as plain properties instead: `sweeps.maxdim = 30` broadcasts (matching
mpscpp3/chain_session.h's `sweeps.maxdim() = maxm_;`), and
`sweeps.maxdim[i]` (0-based) reads one sweep's value.
"""


class Sweeps:
    def __init__(self, nsweep):
        self.nsweep = nsweep
        self._maxdim = [1] * nsweep
        self._cutoff = [0.0] * nsweep
        self._noise = [0.0] * nsweep
        self._niter = [2] * nsweep

    def _prop(name):
        def getter(self):
            return getattr(self, "_" + name)

        def setter(self, value):
            if isinstance(value, (list, tuple)):
                if len(value) != self.nsweep:
                    raise ValueError("Sweeps.{}: expected {} values".format(name, self.nsweep))
                setattr(self, "_" + name, list(value))
            else:
                setattr(self, "_" + name, [value] * self.nsweep)
        return property(getter, setter)

    maxdim = _prop("maxdim")
    cutoff = _prop("cutoff")
    noise = _prop("noise")
    niter = _prop("niter")
    del _prop

    def setnoise(self, i, value):
        """0-based sweep index, matching chain_session.h's
        `for (int i=nsweeps_/2;i<nsweeps_;i++) sweeps.setnoise(i,0.0);`"""
        self._noise[i] = value

    def at(self, i):
        """0-based: (maxdim, cutoff, noise, niter) for sweep i."""
        return self._maxdim[i], self._cutoff[i], self._noise[i], self._niter[i]
