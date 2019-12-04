#!/bin/sh
#

 echo "simulation of NEB calculation: sh run_neb.sh"

date

sh ../../../environment_variables

mpijob -np 432  neb.x -npool 18   -inp  H3S_NEB.in   >  H3S_NEB.out 

echo  
echo "##### done: run.sh  #####"

date

