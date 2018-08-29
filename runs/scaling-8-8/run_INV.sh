#!/bin/bash

# starting Inversion, assuming FW has already been run.
mpirun -np 64 ./IFOS_INV ./INV_8_8.json
