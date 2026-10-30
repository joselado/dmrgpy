import numpy as np, time
n, m = 102400, 34
rng = np.random.default_rng(0)
Bc = np.empty((n, m), complex); Br = np.empty((m, n), complex)
X = rng.standard_normal((m, n)) + 1j*rng.standard_normal((m, n))
Bc[:] = X.T; Br[:] = X
w0 = rng.standard_normal(n) + 1j*rng.standard_normal(n)
def col(k):
    Qk = Bc[:, :k]; w = w0
    h = Qk.conj().T @ w; w = w - Qk @ h; h2 = Qk.conj().T @ w; w = w - Qk @ h2; return np.linalg.norm(w)
def row(k):
    Qk = Br[:k]; w = w0
    h = np.conj(Qk @ np.conj(w)); w = w - h @ Qk; h2 = np.conj(Qk @ np.conj(w)); w = w - h2 @ Qk; return np.linalg.norm(w)
for f in (col, row):
    t = time.time()
    for k in range(1, m + 1): f(k)
    print(f.__name__, "%.1f ms per iteration (mean over k=1..%d)" % ((time.time() - t) / m * 1e3, m))
