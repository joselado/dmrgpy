import sys; sys.path.insert(0, sys.path[0])
from cases import *
from dmrgpy.pyitensor import clik
from dmrgpy.pyitensor.mpsalgebra import randomMPS, inner
bc, lc = bh(5, 4)
S = lc.doubled._session; G = S._mpo(lc.L.to_terms())
rng = np.random.default_rng(3)
states = [randomMPS(S.sites, m, rng=rng) for m in (64, 20, 64, 5, 64, 64, 33)]
km = clik.KrylovMatrices(G, False)
for s in states:
    s.normalize(); km.add(s)
    print("chunks", len(km._chunks(s)), flush=True)
M = np.array([[inner(a, b) for b in states] for a in states])
Gm = np.array([[inner(a, G, b) for b in states] for a in states])
print("M err %.1e  G err %.1e (scale %.1e)" % (np.max(np.abs(M - km.M)), np.max(np.abs(Gm - km.Gm)), np.max(np.abs(Gm))))
