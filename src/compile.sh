#!/bin/bash
# Compile all modules and programs so that they're ready to run.

start=$(date +%s)
echo "Compiling modules..."
gfortran -O3 -J bin/ -o bin/kind_parameters.o -c src/kind_parameters.f90

echo "Compiling gillespie.f90"
gfortran -O3 -o bin/gillespie.o -c src/gillespie.f90 -I bin/
gfortran -O3 -o bin/gillespie.out bin/gillespie.o bin/kind_parameters.o
end=$(date +%s)
echo "Elapsed Time: $(($end-$start)) seconds."
