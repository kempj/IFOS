#!/bin/bash

# starting IFOS for forward modeling
mpirun -np 32 ./IFOS_FW ./FW_8_4.json



#rm jacobian/toy_example/*toy_example*.*.*.*
#rm model/*.bin.*.*
#rm model/toy_example/*.bin.*.*
#rm *.bin.*.*
