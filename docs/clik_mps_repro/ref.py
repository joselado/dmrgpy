import numpy as np, scipy.linalg as sl
sx=np.array([[0,1],[1,0]])/2; sy=np.array([[0,-1j],[1j,0]])/2; sz=np.diag([.5,-.5]); I=np.eye(2)
def op(o,i):
    m=np.array([[1.]])
    for k in range(4): m=np.kron(m,o if k==i else I)
    return m
S=lambda o,i: op(o,i)
h0=0.3*S(sz,0)-0.2*S(sz,3)+0.4*S(sx,1)+0.5*S(sz,2)+0.1*S(sx,2)
w,v=np.linalg.eigh(h0); psi0=v[:,0]
h=sum(S(a,i)@S(a,i+1) for i in range(3) for a in (sx,sy,sz))+0.5*S(sx,0)+0.7j*S(sz,0)
O=S(sz,2); U=sl.expm(-1j*0.05*h)
psi=psi0.copy(); un=[]; nr=[]
for t in range(40):
    un.append(np.vdot(psi,O@psi)); nr.append(np.vdot(psi,O@psi)/np.vdot(psi,psi))
    psi=U@psi
print("unnormalized",np.round(np.real(un[::8]),4)); print("normalized  ",np.round(np.real(nr[::8]),4))
