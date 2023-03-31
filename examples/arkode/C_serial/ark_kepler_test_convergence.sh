#!/bin/bash

# generate reference solution - use 8th order ERK method with tiny timestep
./ark_kepler 1 0 8 0.000001

orders=(1 2 22 222 3 33 4 44 5 6 8 10)
dts=(0.5 0.1 0.01 0.001 0.0001)
for order in ${orders[@]}; 
do 
  for dt in ${dts[@]};
  do 
    # ./ark_kepler 0 0 $order $dt 0 # no compensated sums
    ./ark_kepler 0 0 $order $dt 1 # compensated sums
  done
done

# plot
./ark_kepler_plot_order_work.py
