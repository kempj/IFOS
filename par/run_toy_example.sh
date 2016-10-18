#!/bin/bash

###############################################################
#  running forward modeling for calculation of observed data  #
###############################################################

# compiling all libraries and  IFOS
make install MODEL=../genmod/toy_example_true.c

#-------------------------------------------------------------------------------------
# Insert the trench "Ettlinger Linie" into the toy_example model 
# in the file ../genmod/toy_example_true_trench.c the size and position can be adjusted
# Either an elliptic trench or a v-shaped trench can be choosen

#make install MODEL=../genmod/toy_example_true_trench.c
#---------------------------------------------------------------------------------------

# starting IFOS for forward modeling
# (depending on the MPI package you maybe have to adjust the following programme call, 
#  e.g. when using openMPI you do not have to use the 'lamboot' command)
#lamboot
mpirun -np 4 nice -19 ../bin/IFOS2D in_and_out/toy_example/toy_example_FW.json | tee in_and_out/toy_example/toy_example_FW.out

###############################################################
#                    running the inversion                    #
###############################################################

# compiling IFOS
make install MODEL=../genmod/toy_example_start.c

# starting IFOS
#lamboot
mpirun -np 4 nice -19 ../bin/IFOS2D in_and_out/toy_example/toy_example_INV.json | tee in_and_out/toy_example/toy_example_INV.out
