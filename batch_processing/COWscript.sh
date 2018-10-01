#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem=512
#SBATCH --time=10:00:00

# Load relevant modules
module load intel parallel

# Run the executable
./run_parallel.sh

# Indicate normal completion
exit 0
