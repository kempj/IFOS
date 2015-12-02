#!/bin/bash

###############################################################
#  running forward modeling for calculation of observed data  #
###############################################################

# compiling all libraries and  DENISE
# make clean
make denise MODEL=../genmod/toy_example_true.c
pause
# starting DENISE for forward modeling
# (depending on the MPI package you maybe have to adjust the following programme call, 
#  e.g. when using openMPI you do not have to use the 'lamboot' command)
#lamboot
mpirun -np 4 ../bin/denise in_and_out/toy_example/toy_example_FW.json | tee in_and_out/toy_example/toy_example_FW.out

# the forward modeled data have to be renamed for the inversion
for (( i=1; i <= 5; i++ )) ; do

     mv su/measured_data/toy_example_vx.su.shot${i}.it1 su/measured_data/toy_example_x.su.shot${i}
     mv su/measured_data/toy_example_vy.su.shot${i}.it1 su/measured_data/toy_example_y.su.shot${i}
     mv su/measured_data/toy_example_vz.su.shot${i}.it1 su/measured_data/toy_example_z.su.shot${i}

done

exit

###############################################################
#                    running the inversion                    #
###############################################################

# compiling DENISE
make clean
make denise MODEL=../genmod/toy_example_start.c

# starting DENISE
#lamboot
mpirun -np 4 ../bin/denise in_and_out/toy_example/toy_example_INV.json | tee in_and_out/toy_example/toy_example_INV.out
