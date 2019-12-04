#!/bin/sh
#

 echo " run the second NEB calcultion with larger cutoff:  sh  gen_second_neb.sh  "

###  awk ' {print $0, "   0 ", "   0 ", "   0 " } '  H3S_NEB.xyz  > XYZ.dat

cp H3S_NEB.crd  XYZ.dat

f="H3S_NEB.in"

fin="H3S_initial_2.in"

for ((i=1; i <= 18 ; i++))
do
        d=`expr 18 \* $i`

        if [ $i -eq 1 ]; then
                head -n 50  $fin  >  $f
        elif [ $i -lt 18 ] ; then
                echo 'INTERMEDIATE_IMAGE' >> $f
        else
                echo 'LAST_IMAGE' >> $f
        fi

        echo  'ATOMIC_POSITIONS {crystal}' >> $f

        head -n $d XYZ.dat | tail -n 16  >>  $f

        echo $i
done

echo 'END_POSITIONS' >> $f 

grep -A1 'K_POINTS' $fin >> $f

grep -A3 'CELL_PARAMETERS' $fin >> $f

echo 'END_ENGINE_INPUT' >> $f
echo 'END' >> $f


# tail -n 10 H3S.in  >> $f

echo "##### done #####"


