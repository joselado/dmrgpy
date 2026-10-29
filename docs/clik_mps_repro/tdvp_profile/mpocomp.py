from common import *
def mpo_arrays(G):
    n=G.length(); out=[]
    for i in range(1,n+1):
        T=G.A(i); link=[ind for ind in T.inds if ind.hastags("Link")]
        site0=[ind for ind in T.inds if ind.hastags("Site") and ind.plev==0][0]
        site1=[ind for ind in T.inds if ind.hastags("Site") and ind.plev==1][0]
        order=([link[0]] if i>1 else [])+[site0,site1]+([link[-1]] if i<n else [])
        a=T.transpose_to(order)
        if i==1: a=a[None]
        if i==n: a=a[...,None]
        out.append(np.array(a))
    return out
def compress(Ws, cutoff=1e-13):
    Ws=[w.copy() for w in Ws]; n=len(Ws)
    for i in range(n-1):
        a=Ws[i]; sh=a.shape; q,r=np.linalg.qr(a.reshape(-1,sh[-1]))
        Ws[i]=q.reshape(sh[:-1]+(q.shape[1],)); Ws[i+1]=np.tensordot(r,Ws[i+1],axes=(1,0))
    for i in range(n-1,0,-1):
        a=Ws[i]; sh=a.shape; u,s,vh=np.linalg.svd(a.reshape(sh[0],-1),full_matrices=False)
        k=int(np.sum(s>cutoff*s[0])); u,s,vh=u[:,:k],s[:k],vh[:k]
        Ws[i]=vh.reshape((k,)+sh[1:]); Ws[i-1]=np.tensordot(Ws[i-1],u*s,axes=(Ws[i-1].ndim-1,0))
    return Ws
def mpo_inner(A,B):
    E=np.ones((1,1),complex)
    for a,b in zip(A,B):
        E=np.einsum('xy,xsta,ystb->ab',E,a.conj(),b)
    return E[0,0]
if __name__=="__main__":
    S, G, psi = lindblad_case()
    Ws=mpo_arrays(G)
    C=compress(Ws)
    print("orig",[w.shape[0] for w in Ws]); print("comp",[w.shape[0] for w in C])
    d=mpo_inner(Ws,Ws)+mpo_inner(C,C)-2*mpo_inner(Ws,C)
    print("rel diff", abs(d)/abs(mpo_inner(Ws,Ws)))
