# Bare fermion ladder operators for the Jordan-Wigner-dressed MPO route of
# mpo.py::get_MPO. A parity-odd operator (a single C or Cdag, or any term
# with an odd number of them) cannot go through ITensorMPS's OpSum with the
# raw names: opsum_to_mpo errors "Parity-odd fermionic terms not yet
# supported". For those mpo.py sends the form multioperator.to_terms()
# gives instead, C_i = F_1...F_{i-1} A_i, where A/Adag are plain local
# matrices with NO fermion string, so OpSum treats the term as bosonic and
# the explicit F factors carry the sign. This is the form v3 and "python"
# consume.
#
# Only SiteType"Fermion" needs them. SiteType"Electron" already defines
# Aup/Adagup/Adn/Adagdn (no fermion string), F, Fup and Fdn in ITensors
# itself, so they are deliberately not redefined here.

ITensors.op(::OpName"A", ::SiteType"Fermion") = [0.0 1.0; 0.0 0.0]
ITensors.op(::OpName"Adag", ::SiteType"Fermion") = [0.0 0.0; 1.0 0.0]
ITensors.has_fermion_string(::OpName"A", ::SiteType"Fermion") = false
ITensors.has_fermion_string(::OpName"Adag", ::SiteType"Fermion") = false
