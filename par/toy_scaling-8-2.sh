#!/bin/bash

###############################################################
#  running forward modeling for calculation of observed data  #
###############################################################

rm model_true/mod_toy_example*
rm model/mod_toy_example*
rm model/toy_example/mod_toy_example*
rm jacobian/toy_example/jac_toy_example*
rm su/toy_example/toy_example*
rm su/measured_data/toy_example*

# compiling all libraries and  IFOS
make install MODEL_AC=../genmod/toy_example_ac_true.c

# starting IFOS for forward modeling
mpirun -np 16 nice -19 ../bin/IFOS2D in_and_out/toy_scaling/FW_8_2.json

###############################################################
#                    running the inversion                    #
###############################################################

# compiling IFOS
make install MODEL_AC=../genmod/toy_example_ac_start.c

# starting IFOS
mpirun -np 16 nice -19 ../bin/IFOS2D in_and_out/toy_scaling/INV_8_2.json

make clean

rm jacobian/toy_example/*toy_example*.*.*.*
rm model/*.bin.*.*
rm model/toy_example/*.bin.*.*
rm *.bin.*.*
