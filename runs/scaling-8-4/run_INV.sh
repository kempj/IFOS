#!/bin/bash

# starting Inversion, assuming FW has already been run.
mpirun -np 32 ./IFOS_INV ./INV_8_4.json

