#!/bin/bash

#SBATCH -A snic2019-1-9
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 12:00:00

./results.sh > timings.out
