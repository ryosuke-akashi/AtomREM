#!/bin/sh
#

###  awk ' {print $0, "   0 ", "   0 ", "   0 " } '  H3S_NEB.xyz  > XYZ.dat

cp H3S_NEB.xyz  XYZ.dat

f="H3S_test.in"

for ((i=1; i <= 18 ; i++))
do
        d=`expr 18 \* $i`

        if [ $i -eq 1 ]; then
                head -n 50 H3S.in > $f
        elif [ $i -lt 18 ] ; then
                echo 'INTERMEDIATE_IMAGE' >> $f
        else
                echo 'LAST_IMAGE' >> $f
        fi

        echo  'ATOMIC_POSITIONS {angstrom}' >> $f

        head -n $d XYZ.dat | tail -n 16  >>  $f

        echo $i
done

tail -n 10 H3S.in  >> $f

echo "##### done #####"


