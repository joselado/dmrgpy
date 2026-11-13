# New methods

This folder holds methods that are new in dmrgpy, in the sense that
they are not a direct implementation of one published paper. Each one
has its own folder with the derivation, the papers it rests on, the
checks it passes with their measured numbers, and what it does not
contain, so that the method can be checked again later. The user guide
describes a method only once its checks pass and it is exposed in the
public API.

- `kondo_leading_log/`: the leading-log resummation of the third-order
  STM Kondo spectrum for an impurity with any spectrum, at any
  temperature, in `kondospectrumtk/leadinglog.py` and
  `kondospectrumtk/thermalkernel.py`.
- `kondo_krylov_thermal/`: the thermal third-order STM Kondo spectrum
  without full diagonalization, from one block Krylov space per occupied
  state, in `kondospectrumtk/krylovkondo.py`
  (`get_kondo_spectrum(mode="ED", submode="KRYLOV")`).
- `kondo_nrg_mps/`: the STM spectrum of one spin or a chain of spins on
  a metal by NRG with the complete basis and the full density matrix,
  from the exact first shell or from the bond basis of a DMRG head, in
  `nrgspectral.py` and `nrgmps.py`
  (`ysr.get_stm_spectrum(submode="NRG")`).
