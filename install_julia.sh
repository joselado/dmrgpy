#!/bin/sh

# install Itensor julia
julia install/juliainstall.jl # execute the julia installation
# install python dependences
pip install julia
