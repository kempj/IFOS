#!/bin/bash

# starting Inversion, assuming FW has already been run.
mpirun -np 12 ./IFOS_INV ./INV_1_12.json

