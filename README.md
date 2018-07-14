## DMRGPY ##

# Summary #

This library is a Python wrapper to compute spin chains using density
matrix renormalization group as implemented in ITensor.

Several examples can be found in the examples folder.

# Disclaimer #

This library is still under heavy development

# How to install #

The script install.sh will compile both ITensor and a C++ program
that uses it it. Afterwards, it is only required to add to the .bashrc
the following line

export DMRGROOT=PATH_TO_DMRGPY"/src"

# Capbilities #
- Ground state energies
- Gap
- Excited states
- Correlation functions
- Dynamical correlation functions with Kernel polynomial method
