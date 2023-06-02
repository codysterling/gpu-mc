#!/usr/bin/env sh

# Running the serial code, arguments are: natoms nsteps wsteps
for i in 125 1000 3375 8000
do
  echo "Now running with $i atoms"
  for j in {1..5}
  do
    echo "Doing run $j:"
    OMP_NUM_THREADS=$1 ./cody-mc-cpu-parallel.out $i 1000 10
    OMP_NUM_THREADS=$1 ./cody-mc-cpu-parallel.out $i 100000 1000
    OMP_NUM_THREADS=$1 ./cody-mc-cpu-parallel.out $i 500000 10000
    OMP_NUM_THREADS=$1 ./cody-mc-cpu-parallel.out $i 1000000 10000
    if [ $j -lt 5 ]; then echo ""; fi
  done
  echo "Finished run of $i atoms"
  echo ""
done
