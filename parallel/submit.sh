#!/bin/bash

#SBATCH -A snic2019-1-9
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 12:00:00

for i in 1 2 4 8 16 32
do
  ./results.sh $i > timings-${i}.out
done
