# compile Itensor library
cd src/dmrgpy/mpscpp/ITensor-master
cp options.save options.mk
make clean 
make
cd -

# compile DMRG program
cd src/dmrgpy/mpscpp
make clean 
make
mv mpscpp mpscpp.x
cd -

# compile fortran routines of classic dmrg
cd src/dmrgpy/pychain/fortran
bash compile_fortran.sh
cd -


python addsystem.py # add route to .bashrc
