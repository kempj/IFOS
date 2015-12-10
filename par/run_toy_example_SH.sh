#!/bin/bash

###############################################################
#  running forward modeling for calculation of observed data  #
###############################################################

# compiling all libraries and  DENISE
make clean
make denise MODEL=../genmod/toy_example_true.c MODEL_EL=../genmod/toy_example_elastic_true.c

# starting DENISE for forward modeling
# (depending on the MPI package you maybe have to adjust the following programme call, 
#  e.g. when using openMPI you do not have to use the 'lamboot' command)
#lamboot
mpirun -np 4 ../bin/denise in_and_out/toy_example/toy_example_FW_SH.json | tee in_and_out/toy_example/toy_example_FW_SH.out

# the forward modeled data have to be renamed for the inversion
for (( i=1; i <= 5; i++ )) ; do
     mv su/measured_data/toy_example_vz.su.shot${i}.it1 su/measured_data/toy_example_z.su.shot${i}

done


###############################################################
#                    running the inversion                    #
###############################################################

# compiling DENISE
make clean
make denise MODEL=../genmod/toy_example_start.c MODEL_EL=../genmod/toy_example_elastic_start.c

# starting DENISE
#lamboot
mpirun -np 4 ../bin/denise in_and_out/toy_example/toy_example_INV_SH.json | tee in_and_out/toy_example/toy_example_INV_SH.out
