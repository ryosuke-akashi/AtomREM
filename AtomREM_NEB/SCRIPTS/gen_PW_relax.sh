#!/bin/sh
#
 echo 
 echo "run generation of PW relaxation for cell parameters after NEB: sh  gen_PW_relax.sh"
 
 awk ' {print $0, "   0 ", "   0 ", "   0 " } '  H3S_NEB.crd  > XYZ.dat

fin='H3S.in'


for ((i=1; i <= 18 ; i++))
do

echo "Image number  " $i

d=`expr 18 \* $i`

echo $d

dir="IM_"$i

mkdir ${dir}

f=${dir}"/H3S.in"


head -n 37 $fin  >  $f
head -n $d XYZ.dat | tail -n 16  >>  $f

grep -A1 'K_POINTS' $fin >> $f

grep -A3 'CELL_PARAMETERS' $fin >> $f

echo "" >> $f 

sed "s/IM_1/IM_$i/g" $f > 1.txt 

mv 1.txt $f  

## cp run.sh ${dir}

### sh run.sh 

## sbatch -n24 -p batch submit.sh

echo "Image number  " $i  "  is done. "
echo 

done

echo "##### done #####"


