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

python addsystem.py # add route to .bashrc
