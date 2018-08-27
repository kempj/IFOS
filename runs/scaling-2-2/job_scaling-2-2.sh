#!/bin/bash
#SBATCH -J IFOS2D-2-2      # job name
#SBATCH -o IFOS2D-2-2.o%j  # output and error file name (%j expands to jobID)
#SBATCH -N 2                    # total number of nodes requested
#SBATCH --ntasks-per-node 2    # total number of processors on each node
#SBATCH --mem 16000              # total memory request in MB
#SBATCH -p compute                  # queue (partition)
#SBATCH -t 30:00:00             # run time (hh:mm:ss)
#SBATCH --mail-user=jekemp@pvamu.edu  # change to your email address
#SBATCH --mail-type=begin       # email me when the job starts
#SBATCH --mail-type=end         # email me when the job finishes

./toy_scaling-2-2.sh

