#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem=512
#SBATCH --time=10:00:00

# Load relevant modules
source ~/.bashrc

# Run the executable
./run.sh

# Indicate normal completion
exit 0
