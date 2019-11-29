#!/bin/sh
#

 awk ' {print $0, "   0 ", "   0 ", "   0 " } '  H3S_NEB.xyz  > XYZ.dat

for ((i=1; i <= 18 ; i++))
do

d=`expr 18 \* $i`

echo $d

dir="IM_"$i

mkdir ${dir}

f=${dir}"/H3S.in"


head -n 37 H3S.in > $f
head -n $d XYZ.dat | tail -n 16  >>  $f
tail -n 7 H3S.in  >> $f

cp submit.sh ${dir}

cd ${dir}



sbatch -n24 -p batch submit.sh

cd ../


### qsub stat_$i.sh 

## f="PATHS/path_"$i".dat"

## make a simulation 

### cp path.dat $f

echo $i

done

echo "##### done #####"


