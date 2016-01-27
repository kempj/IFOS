#!/bin/bash

###############################################################
#  running forward modeling for calculation of observed data  #
###############################################################

# compiling all libraries and  IFOS
make clean
make IFOS2D MODEL=../genmod/toy_example_true.c MODEL_EL=../genmod/toy_example_elastic_true.c

# starting IFOS for forward modeling
# (depending on the MPI package you maybe have to adjust the following programme call, 
#  e.g. when using openMPI you do not have to use the 'lamboot' command)
#lamboot
mpirun -np 4 ../bin/IFOS2D in_and_out/toy_example/toy_example_FW_SH.json | tee in_and_out/toy_example/toy_example_FW_SH.out


###############################################################
#                    running the inversion                    #
###############################################################

# compiling IFOS
make clean
make IFOS2D MODEL=../genmod/toy_example_start.c MODEL_EL=../genmod/toy_example_elastic_start.c

# starting IFOS
#lamboot
mpirun -np 4 ../bin/IFOS2D in_and_out/toy_example/toy_example_INV_SH.json | tee in_and_out/toy_example/toy_example_INV_SH.out
