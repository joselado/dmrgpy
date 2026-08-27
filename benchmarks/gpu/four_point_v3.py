"""The same four-point tensor on the compiled ITensor v3 backend, as the
baseline the batched kernel is measured against.

Separate script, and a separate process, because it needs a checkout whose
`mpscpp3` extension is actually compiled -- which is not the pure-Python
source snapshot the batched runs use.
"""
import argparse, json, sys, time


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True)
    ap.add_argument("--n", type=int, default=30)
    ap.add_argument("--maxm", type=int, nargs="+", default=[20, 40])
    ap.add_argument("--nsweeps", type=int, default=4)
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    sys.path.insert(0, a.src)
    import numpy as np
    from dmrgpy import fermionchain

    rows = []
    for maxm in a.maxm:
        fc = fermionchain.Fermionic_Chain(a.n, itensor_version=3)
        h = 0
        for i in range(a.n - 1):
            h = h + fc.Cdag[i] * fc.C[i + 1] + fc.Cdag[i + 1] * fc.C[i]
            h = h + 1.0 * fc.N[i] * fc.N[i + 1]
        fc.set_hamiltonian(h)
        fc.maxm = maxm
        fc.nsweeps = a.nsweeps
        t0 = time.time()
        e0 = fc.gs_energy()
        wf = fc.get_gs()
        tgs = time.time() - t0
        t0 = time.time()
        ct = wf.get_four_correlation_tensor()
        t = time.time() - t0
        print("n=%d maxm=%d E0=%.10f gs=%.1fs  v3 four-point %8.2f s  checksum %.10f"
              % (a.n, maxm, np.real(e0), tgs, t, float(np.sum(np.abs(ct)))), flush=True)
        rows.append({"n": a.n, "maxm": maxm, "e0": float(np.real(e0)),
                     "gs_time": tgs, "v3": t,
                     "checksum": float(np.sum(np.abs(ct)))})
    if a.json:
        with open(a.json, "w") as f:
            json.dump(rows, f, indent=1)


if __name__ == "__main__":
    main()
