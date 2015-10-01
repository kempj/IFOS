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

# compiling all libraries and  DENISE
make clean
make denise MODEL_AC=../genmod/toy_example_ac_true.c

# starting DENISE for forward modeling
mpirun -np 4 nice -19 ../bin/denise in_and_out/toy_example/toy_example_ac_FW.json | tee in_and_out/toy_example/toy_example_ac_FW.out
# mpirun -np 4 xterm -e gdb --args ../bin/denise in_and_out/toy_example/toy_example_ac_FW.json | tee in_and_out/toy_example/toy_example_ac_FW.out
# the forward modeled data have to be renamed for the inversion
for (( i=1; i <= 5; i++ )) ; do

     mv su/measured_data/toy_example_ac_p.su.shot${i}.it1 su/measured_data/toy_example_ac_p.su.shot${i}
     mv su/measured_data/toy_example_ac_vx.su.shot${i}.it1 su/measured_data/toy_example_ac_x.su.shot${i}
     mv su/measured_data/toy_example_ac_vy.su.shot${i}.it1 su/measured_data/toy_example_ac_y.su.shot${i}

done

###############################################################
#                    running the inversion                    #
###############################################################

# compiling DENISE
make clean
make denise MODEL_AC=../genmod/toy_example_ac_start.c

# starting DENISE
mpirun -np 4 nice -19 ../bin/denise in_and_out/toy_example/toy_example_ac_INV.json | tee in_and_out/toy_example/toy_example_ac_INV.out

make clean

# rm jacobian/toy_example/*.old.*.*
# rm model/toy_example/*.bin.*.*
