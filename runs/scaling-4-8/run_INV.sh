#!/bin/bash


# starting Inversion, assuming FW has already been run.
mpirun -np 16 ./IFOS_INV ./INV_4_8.json
