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
