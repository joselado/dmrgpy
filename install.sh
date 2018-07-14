# compile Itensor library
cd src/mpscpp/ITensor-master
cp options.save options.mk
make clean 
make
cd -

# compile DMRG program
cd src/mpscpp
make clean 
make
cd -

