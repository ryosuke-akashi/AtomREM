#!/bin/sh

## tail -n 7 STRUCTURES/atoms_0.dat



##for ((j=1; j < 41; j++))

for j in 7 12 14 23 27 28 31 32 33 34 35 36 37 38 40
do
let i=j-1

f="STRUCTURES/atoms_"$i".dat"
d=$j"th_ROBBER/"

mkdir $d

f2=$d"/atoms.dat"

##echo $f
#echo $f2

tail -n 7 $f > $f2

echo $j
echo $d



cp param_3N.in $d
cp submit.sh $d
cp a.out  $d


##echo $f

cd $d 
qsub submit.sh 
cd ../


done








