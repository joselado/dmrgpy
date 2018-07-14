## DMRGPY ##

# Summary #

This is a Python library to compute spin chains using density
matrix renormalization group as implemented in ITensor. Most
of the camputations can be performed both with DMRG and exact
diagonalization for small systems, which allows to benchmark the
results.

Several examples can be found in the examples folder.

# Disclaimer #

This library is still under heavy development

# How to install #

The script install.sh will compile both ITensor and a C++ program
that uses it. Afterwards, it is only required to add to the .bashrc
the following line

export DMRGROOT=PATH_TO_DMRGPY"/src"

The whole Python library will use that environmental variable to
find the necesary programs.

# Capabilities #
- Ground state energy
- Gap
- Excited states
- Correlation functions
- Dynamical correlation functions with Kernel polynomial method
