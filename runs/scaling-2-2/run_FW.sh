#!/bin/bash

###############################################################
#  running forward modeling for calculation of observed data  #
###############################################################

#rm model_true/mod_toy_example*
#rm model/mod_toy_example*
#rm model/toy_example/mod_toy_example*
#rm jacobian/toy_example/jac_toy_example*
#rm su/toy_example/toy_example*
#rm su/measured_data/toy_example*

# starting IFOS for forward modeling
mpirun -np 4 nice -19 ./IFOS_FW ./FW_2_2.json



#rm jacobian/toy_example/*toy_example*.*.*.*
#rm model/*.bin.*.*
#rm model/toy_example/*.bin.*.*
#rm *.bin.*.*
