from common import *
from mpocomp import mpo_arrays
S, G, psi = lindblad_case()
Ws=mpo_arrays(G)
W=Ws[3]
np.set_printoptions(linewidth=250)
# classify each (a,b) by its operator matrix (normalized)
D1,_,_,D2=W.shape
pat=np.zeros((D1,D2),dtype=object)
ops=[]
for a in range(D1):
    for b in range(D2):
        M=W[a,:,:,b]
        if np.abs(M).max()<1e-14: pat[a,b]="."; continue
        # identify by support signature
        sig=tuple(map(tuple,np.argwhere(np.abs(M)>1e-14)))
        if sig not in ops: ops.append(sig)
        pat[a,b]=str(ops.index(sig))
for a in range(D1): print(" ".join("%2s"%x for x in pat[a]))
print(len(ops),"distinct support patterns")
for k,s in enumerate(ops): print(k,len(s),s[:6])
