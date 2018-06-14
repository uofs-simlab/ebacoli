#!/usr/bin/env bash

# script to perform multiple C1 value runs

driverTemplate=TEMPLATE_driver-CDM-frameoutput.f95
echo "############################################################################"
echo "############################################################################"
echo "Generating ${PWD##*/} data"

C1_vals="1.0D0 1.1D0"

for C1 in ${C1_vals}
do
    echo "############################################################################"
    echo "Running C1=${C1}"
    echo "############################################################################"
    # replace C1 value
    sed  "s/#C1#/${C1}/" < $driverTemplate > ${driverTemplate#TEMPLATE_}
    make run-all

    dirName=C1_${C1}
    mkdir -p ${dirName}
    mv Bsplines000* ${dirName}
    mv heatmap*.eps ${dirName}
done
