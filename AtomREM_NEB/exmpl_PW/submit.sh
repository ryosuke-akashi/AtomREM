#!/bin/sh
#QSUB -queue i9acc 
#QSUB -node 9
#QSUB -mpi  216
#PBS -l walltime=00:30:00
#

date
cd ${PBS_O_WORKDIR}

. ../../../environment_variables

mpijob -np 216  pw.x -npool 9  < H3S.in   >  H3S.out 

echo  
echo "##### done #####"

date

