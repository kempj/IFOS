#!/bin/bash

#rm model/mod_toy_example*
#rm model/toy_example/mod_toy_example*
#rm jacobian/toy_example/jac_toy_example*
#rm su/toy_example/toy_example*
#rm su/measured_data/toy_example*

# starting Inversion, assuming FW has already been run.
mpirun -np 16 ./IFOS_INV ./INV_8_2.json


#rm jacobian/toy_example/*toy_example*.*.*.*
#rm model/*.bin.*.*
#rm model/toy_example/*.bin.*.*
#rm *.bin.*.*
