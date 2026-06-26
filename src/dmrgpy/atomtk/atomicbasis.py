import numpy as np


def get_r2h_d(orbs):
    """"Unitary transformation from d orbitals to spherical harmonics"""
    # we assume input cartesian basis dz2,dxz,dyz,dx2y2,dxy
    # spherical basis is -2,-1,0,1,2
    T1 = np.array([
    [0,   0,              0,              1/np.sqrt(2),  -1j/np.sqrt(2)],
    [0,  -1/np.sqrt(2),  +1j/np.sqrt(2),       0,              0],
    [1,   0,              0,                   0,              0],
    [0,  1/np.sqrt(2),  +1j/np.sqrt(2),        0,              0],
    [0,   0,                0,            1/np.sqrt(2),   +1j/np.sqrt(2)]
    ])
    # now the mapping between the orbitals
    odict = dict()
    odict["dz2"] = 0 # orbital
    odict["dxz"] = 1 # orbital
    odict["dyz"] = 2 # orbital
    odict["dx2y2"] = 3 # orbital
    odict["dxy"] = 4 # orbital
    # now the permutation
    T2 = np.zeros((5,5))
    for io in range(len(orbs)):
        ii = orbs[io] # input orbital
        jj = odict[ii] # output index
#        print(ii,io,jj)
        T2[io,jj] = 1.0 # permutation
#    print(T2)
#    print(T2.T@np.array([0.,0.,1.,0.,0.]))
    return T1@T2.T # output



def get_T(orbs):
    return get_r2h_d(orbs)




from .. import fermionchain


def generate_atom(orbs=None, # orbitals
                  tij=None,  # single particle spinless hopping
                          U=4., # interaction 
                          B=[0.,0.,0.], # magnetic field
                          soc=0., # SOC
                          J=0.5, # Hund's coupling
                          Ne = 5, # number of electrons
                          lamb_Ne = 200, # lagrange multiplier
                          ):
    """Build the Hamiltonian"""
    # get transformation from real orbital to spherical harmonics
    if tij is None: raise
    if orbs is None: raise
#    orbs = ["dz2","dxz","dyz","dx2y2","dxy"]
    T = get_T(orbs)

    # angular momenta operators
    Lzp = np.diag([2, 1, 0, -1, -2])

    # angular momenta for L=2
    Lxp = (1./ 2) * np.array([
    [0, 2, 0, 0, 0],
    [2, 0, np.sqrt(6), 0, 0],
    [0, np.sqrt(6), 0, np.sqrt(6), 0],
    [0, 0, np.sqrt(6), 0, 2],
    [0, 0, 0, 2, 0]
    ])
    Lyp = (1./ (2j)) * np.array([
    [0, 2, 0, 0, 0],
    [-2, 0, np.sqrt(6), 0, 0],
    [0, -np.sqrt(6), 0, np.sqrt(6), 0],
    [0, 0, -np.sqrt(6), 0, 2],
    [0, 0, 0, -2, 0]
    ])
    Td = np.conjugate(T.T) # dagger of transformation
    lz = Td@Lzp@T # spherical to cartesian
    lx = Td@Lxp@T # spherical to cartesian
    ly = Td@Lyp@T # spherical to cartesian

    # Pauli matrices
    sigma_z = 0.5*np.array([[1, 0],
                            [0, -1]])

    sigma_x = 0.5*np.array([[0, 1],
                            [1, 0]])

    sigma_y = 0.5*np.array([[0, -1j],
                            [1j, 0]])

    #################################################################

    n = len(orbs)  # Orbitals
    fc = fermionchain.Spinful_Fermionic_Chain(n)  # Fermionic chain

    def one2many(ml=None,ms=None):
        """Return the many body representation of the single particle
        operators ml in orbital basis and ms in spin basis"""
        o = 0 # output
        if ml is None: ml = np.identity(n) # identity matrix
        if ms is None: ms = np.identity(2) # identity matrix
        for io in range(n): # loop over spinful orbitals
            for jo in range(n): # loop over spinful orbitals
                cas = [fc.Cdagup[io],fc.Cdagdn[io]] # first operator
                cbs = [fc.Cup[jo],fc.Cdn[jo]] # second operator
                for si in range(2): # loop over spin
                    for sj in range(2): # loop over spin
                        coeff = ml[io,jo]*ms[si,sj] # orbital times spin
                        if abs(coeff)>1e-6: # if non-zero
                            o = o + coeff*cas[si]*cbs[sj] # add this contribution
        return o # return

    # initialize Hamiltonian
    h = 0

    # constrain the density
    N_tot = sum(fc.N) # total electron operator
    h += lamb_Ne*(N_tot-Ne)*(N_tot-Ne) # force a filling
    h = h + h.get_dagger() # Lagrange multiplier for the total electron number

    # Hopping in the Hamiltonian
    h = h + one2many(ml=tij) # hopping in many-body basis

    # add the SOC
    Hsoc = one2many(ml=lx,ms=sigma_x) # XX
    Hsoc += one2many(ml=ly,ms=sigma_y) # YY
    Hsoc += one2many(ml=lz,ms=sigma_z) # ZZ
    h = h + soc*Hsoc # add SOC contribution

    Lx = one2many(ml=lx) # Lx operator
    Ly = one2many(ml=ly) # Ly operator
    Lz = one2many(ml=lz) # Lz operator
    Sx = one2many(ms=sigma_x) # Sx operator
    Sy = one2many(ms=sigma_y) # Sy operator
    Sz = one2many(ms=sigma_z) # Sz operator


    S2 = Sz*Sz + Sx*Sx + Sy*Sy     # total spin (many-body)
    L2 = Lz*Lz + Lx*Lx + Ly*Ly     # total L (many-body)

    # Add interaction terms to Hamiltonian
    h -= U * S2 + J * L2
    
    # Add the external magnetic field
    h = h + B[0]*(2*Sx+Lx)
    h = h + B[1]*(2*Sy+Ly)
    h = h + B[2]*(2*Sz+Lz)

    # save the operators
    fc.STx = Sx
    fc.STy = Sy
    fc.STz = Sz
    fc.LTx = Lx
    fc.LTy = Ly
    fc.LTz = Lz
    fc.LT2 = L2
    fc.ST2 = S2
    # uncomment if you want to check the commutations
#    from commutation import test_commutations
#    test_commutations(fc,[Lx,Ly,Lz],[Sx,Sy,Sz])
    fc.set_hamiltonian(h) # set the Hamiltonian
    return fc









