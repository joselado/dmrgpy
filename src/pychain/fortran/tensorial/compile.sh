f2py2.7 -llapack -c -m tensorialf90 tensorialf90.f90 
#f2py2.7 tensorialf90.f90 -m tensorialf90-h tensorialf90.pyf
cp tensorialf90.so ../../ # copy to main direcuy
