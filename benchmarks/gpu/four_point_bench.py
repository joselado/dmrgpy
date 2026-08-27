"""Four-point correlator tensor: the batched kernel on host vs device.

<Cdag_i C_j Cdag_k C_l> on a spinless fermionic chain, timed for
`itensor_version="python"` in three configurations:

  * ctmode="sweep"    -- the pre-existing per-tuple implementation
  * ctmode="batched"  -- pyitensor/fourpoint.py on NumPy
  * ctmode="batched"  -- the same kernel with its arrays on a device

The ground state is always solved on NumPy: at these bond dimensions the
device loses on DMRG (docs/gpu_cpu_performance.md), and the four-point
kernel converts the MPS once, in `_arrays_lpr`, so it does not need the
engine that produced the wavefunction to live where it does.

Run `--help` for the knobs. Pin BLAS threads before starting (see CLAUDE.md).
"""
import argparse, json, os, sys, time


def build(fc_mod, n, version, maxm, nsweeps):
    fc = fc_mod.Fermionic_Chain(n, itensor_version=version)
    h = 0
    for i in range(n - 1):
        h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
        h = h + 1.0 * fc.N[i] * fc.N[i + 1]
    fc.set_hamiltonian(h)
    fc.maxm = maxm
    fc.nsweeps = nsweeps
    return fc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=None, help="dmrgpy src/ to put first on sys.path")
    ap.add_argument("--n", type=int, default=30)
    ap.add_argument("--maxm", type=int, nargs="+", default=[20, 40, 80])
    ap.add_argument("--nsweeps", type=int, default=4)
    ap.add_argument("--device", action="store_true", help="also time the JAX backend")
    ap.add_argument("--with-sweep", type=int, default=0,
                    help="also time ctmode='sweep' for maxm <= this (0 = never)")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    if a.src:
        sys.path.insert(0, a.src)
    import numpy as np
    from dmrgpy import fermionchain
    from dmrgpy.pyitensor import backend

    rows = []
    for maxm in a.maxm:
        fc = build(fermionchain, a.n, "python", maxm, a.nsweeps)
        t0 = time.time()
        e0 = fc.gs_energy()
        wf = fc.get_gs()
        tgs = time.time() - t0
        chi = max(wf.cpp_handle.A(i).inds[0].dim for i in range(2, a.n + 1))
        row = {"n": a.n, "maxm": maxm, "e0": float(np.real(e0)),
               "gs_time": tgs, "chi": int(chi)}
        print("n=%d maxm=%d chi=%d E0=%.10f gs=%.1fs" % (a.n, maxm, chi, np.real(e0), tgs),
              flush=True)

        t0 = time.time()
        ref = wf.get_four_correlation_tensor(ctmode="batched")
        row["batched_host"] = time.time() - t0
        print("  batched (numpy)  %8.2f s" % row["batched_host"], flush=True)

        if a.with_sweep and maxm <= a.with_sweep:
            t0 = time.time()
            old = wf.get_four_correlation_tensor(ctmode="sweep")
            row["sweep_host"] = time.time() - t0
            row["sweep_agree"] = float(np.max(np.abs(old - ref)))
            print("  sweep   (numpy)  %8.2f s   agree %.2e"
                  % (row["sweep_host"], row["sweep_agree"]), flush=True)

        if a.device:
            backend.set_backend("jax")
            print("  device: %s" % backend.device_info(), flush=True)
            # Twice, and both numbers reported. XLA compiles per (op,
            # shape) and this kernel mints a fresh (B, chi, chi) leading
            # shape for nearly every transfer -- B is the environment
            # batch, which differs per trie node per site -- so the first
            # call carries thousands of one-time compilations. Reporting
            # only that would read as "the device loses" when what it
            # measured was the compiler. The host needs no repeat: NumPy
            # has no compile step, and at maxm=80 a second pass costs
            # minutes for nothing.
            t0 = time.time()
            got = np.asarray(wf.get_four_correlation_tensor(ctmode="batched"))
            row["batched_device_cold"] = time.time() - t0
            t0 = time.time()
            got = np.asarray(wf.get_four_correlation_tensor(ctmode="batched"))
            row["batched_device_warm"] = time.time() - t0
            row["device_agree"] = float(np.max(np.abs(got - ref)))
            print("  batched (jax)    %8.2f s cold  %8.2f s warm   agree %.2e"
                  % (row["batched_device_cold"], row["batched_device_warm"],
                     row["device_agree"]), flush=True)
            print("  speedup device/host: %.2fx cold, %.2fx warm"
                  % (row["batched_host"] / row["batched_device_cold"],
                     row["batched_host"] / row["batched_device_warm"]), flush=True)
            backend.set_backend("numpy")
        rows.append(row)

    if a.json:
        with open(a.json, "w") as f:
            json.dump(rows, f, indent=1)


if __name__ == "__main__":
    main()
