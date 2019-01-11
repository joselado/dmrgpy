#f2py3.4 -llapack -c -m density_matrixf90 density_matrix.f90
#cp density_matrixf90.cpython-34m.so ../pychain/density_matrixf90_v3.so

f2py -llapack -c -m density_matrixf90 density_matrix.f90
cp density_matrixf90.so ../pychain/density_matrixf90.so
