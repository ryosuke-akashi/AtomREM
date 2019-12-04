#!/bin/sh
#QSUB -queue F18acc
#QSUB -node 18
#QSUB -mpi  432
#PBS -l walltime=23:30:00
#
#PBS -N bulkjob


date
cd ${PBS_O_WORKDIR}

### There are files of reaction path from AtomREM simulation

f='inout_AtomREM/atomic_coord_q.lammpstrj'
f2='inout_AtomREM/coordinates_rank_100020.lammpstrj'


#####################################################################################################
 echo "FIRST: NEB calculation with a short cutoff = 20 Ry" 


dir="NEB_LESS"

mkdir $dir

cp $f $dir 
cp $f2 $dir

### 'H3S_initial.in' is a file with cell parameters and coordinates getting after pw relaxation

cp H3S_initial.in $dir 

cp SCRIPTS/extract_images.sh  $dir 
cp SCRIPTS/make_IN.sh         $dir 
cp SCRIPTS/run_neb.sh         $dir

cd $dir

### Extract images from reaction path:
sh extract_images.sh "../$f"  "../$f2" 
### Output file is 'images.dat'

### Make inout file with 18 images for NEB calculation: 
sh make_IN.sh  ## input file is 'images.dat'
#### Output is 'H3S_NEB.in' for NEB calculations
 
###  sh ../../../environment_variables

### mpijob -np 432  neb.x -npool 18  -inp  H3S_NEB.in   >  H3S_NEB.out 

sh run_neb.sh 


cd ../


#####################################################################################################
 echo "SECOND: NEB calculation with larger cutoff = 80 Ry " 
 echo "The results of first NEB are in the dir = NEB_LESS "


dir2='NEB_FIN'

cp -r $dir $dir2

cp SCRIPTS/gen_second_neb.sh  $dir2   
cp H3S_initial_2.in   $dir2
cp SCRIPTS/run_neb.sh         $dir2

cd $dir2

##### 'H3S_NEB.crd' and 'H3S_initial_2.in' are input files for gen_second_neb.sh
sh gen_second_neb.sh 
### 'H3S_NEB.in' is output file to continue the simulation 

## sh ../../../environment_variables

## mpijob -np 432  neb.x -npool 18  -inp  H3S_NEB.in   >  H3S_NEB.out

sh run_neb.sh

cd ../

#####################################################################################################
 echo " THIRD: PW relaxation after NEB "


dir='PW_NEB'

rm -rf $dir 
mkdir   $dir 

cp H3S.in $dir 
cp $dir2/H3S_NEB.crd  $dir 
cp SCRIPTS/gen_PW_relax.sh    $dir 
cp SCRIPTS/run.sh $dir  

cp SCRIPTS/extract_minimum_energies.sh  $dir 


cd $dir 

sh gen_PW_relax.sh 

sh run.sh 


sh extract_minimum_energies.sh 

cp path_rc_relax.dat   ../


cd ../

echo  
echo "##### done #####"

date




