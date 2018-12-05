#!/bin/bash

# Compile ITensor
cd itensor/ITensor-master
make
cd -


# Compile DMRG program
cd itensor
make
cd -



