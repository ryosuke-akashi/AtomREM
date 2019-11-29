#!/bin/sh
#QSUB -queue i18cpu 
#QSUB -node 18
#QSUB -mpi  432
#PBS -l walltime=00:30:00
#

date
cd ${PBS_O_WORKDIR}

dir="NEB_LESS_2"

mkdir $dir

cp H3S.in $dir 

cd $dir

sh ../../../environment_variables

### mpijob   ./a.out < ./param_3N.in > out


mpijob -np 432  neb.x -npool 18  -nd 24  -inp  H3S.in   >  H3S.out 

cd ../


echo  
echo "##### done #####"

date

