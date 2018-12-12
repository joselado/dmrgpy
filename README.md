## DMRGPY ##

# Summary #

This is a Python library to compute quasi-one-dimensional
spin chains and fermionic systems using matrix product states
with the matrix renormalization group as implemented in ITensor. Most
of the camputations can be performed both with DMRG and exact
diagonalization for small systems, which allows to benchmark the
results.

Several examples can be found in the examples folder.

# Disclaimer #

This library is still under heavy development.

# How to install #

The script install.sh will compile both ITensor and a C++ program
that uses it. Afterwards, it is only required to add to the .bashrc
the following line

export DMRGROOT=PATH_TO_DMRGPY"/src"

After this, you can write in your Python script

import os
import sys
sys.path.append(os.environ["DMRGROOT"])

And import the sublibrary that you want, for example

from dmrgpy import spinchain

# Capabilities #
- Ground state energy
- Excitation gap
- Excited states
- Static correlation functions
- Time evolution and measurements
- Dynamical correlation functions computed with the Kernel polynomial method
- Dynamical correlation functions with time dependent DMRG
