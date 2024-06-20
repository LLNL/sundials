#!/bin/bash

for reactor in 0 1 2
do
  for nls in 0 1 2 3 12 13 22 23
  do
    tname="reactor-${reactor}_nls-${nls}"
    echo $tname
    export SUNLOGGER_DEBUG_FILENAME=${tname}.debug.log
    export SUNLOGGER_INFO_FILENAME=${tname}.debug.log
    ./cvBrusselator ${reactor} ${nls} 0 &> ${tname}.log
    mv ./cvBrusselator_solution.txt ${tname}.solution.txt
    mv ./cvBrusselator_ele.txt ${tname}.ele.txt
    ./plot_cvBrusselator.py ${tname}.solution.txt ${tname}.debug.log ${tname}.ele.txt
    mv cvBrusselator.png ${tname}.png
    mv cvBrusselator_ele.png ${tname}_ele.png
  done
done
