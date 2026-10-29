import sys; sys.path.insert(0,"/u/40/ladovj1/data/Documents/programs/dmrgpy/src")
import numpy as np
from dmrgpy import spinchain
from dmrgpy.lindblad import Lindblad_Chain
sc1=spinchain.Spin_Chain([3,3],itensor_version="python")
H1=sc1.Sx[0]*sc1.Sy[1]+0.3*sc1.Sz[0]*sc1.Sz[0]; sc1.set_hamiltonian(H1)
lc=Lindblad_Chain(sc1,[0.5*(sc1.Sx[0]+1j*sc1.Sy[0]), 0.4*sc1.Sy[1]*sc1.Sz[1]])
A=lc.doubled.get_ED_obj().get_operator(lc.L).toarray()
ed=sc1.get_ED_obj(); dense=lambda o: ed.get_operator(o).toarray()
D=9; I=np.eye(D); Hd=dense(H1); B=-1j*(np.kron(Hd,I)-np.kron(I,Hd.T))
for J in lc.jumps:
    Jd=dense(J); JdJ=Jd.conj().T@Jd; B+=np.kron(Jd,Jd.conj())-0.5*np.kron(JdJ,I)-0.5*np.kron(I,JdJ.T)
for k in range(1,6):
    a=np.trace(np.linalg.matrix_power(A,k)); b=np.trace(np.linalg.matrix_power(B,k)); print(k, abs(a-b)/abs(b))
ev=np.linalg.eigvals(B); print("min eigval gap", np.min([abs(x-y) for i,x in enumerate(ev) for y in ev[i+1:]]))
