#!/bin/bash

###############################################################
#  running forward modeling for calculation of observed data  #
###############################################################

# compiling all libraries and  DENISE
make clean
make denise MODEL=../genmod/toy_example_true.c

# starting DENISE for forward modeling
# (depending on the MPI package you maybe have to adjust the following programme call, 
#  e.g. when using openMPI you do not have to use the 'lamboot' command)
#lamboot
mpirun -np 4 nice -19 ../bin/denise in_and_out/toy_example/toy_example_FW.json | tee in_and_out/toy_example/toy_example_FW.out

###############################################################
#                    running the inversion                    #
###############################################################

# compiling DENISE
make clean
make denise MODEL=../genmod/toy_example_start.c

# starting DENISE
#lamboot
mpirun -np 4 nice -19 ../bin/denise in_and_out/toy_example/toy_example_INV.json | tee in_and_out/toy_example/toy_example_INV.out
