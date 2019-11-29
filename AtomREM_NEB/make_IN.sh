#!/bin/sh
#

###  awk ' {print $0, "   0 ", "   0 ", "   0 " } '  H3S_NEB.xyz  > XYZ.dat

cp images.dat XYZ.dat

# cp $1  XYZ.dat

f_in="H3S_initial.in"

f="H3S_NEB.in"

n_atoms=16
n_lines=`expr $n_atoms + 9`

rm $f

for ((i=1; i <= 18 ; i++))
do
        d=`expr $n_lines \* $i`

        if [ $i -eq 1 ]; then
                head -n 50 $f_in > $f
        elif [ $i -lt 18 ] ; then
                echo 'INTERMEDIATE_IMAGE' >> $f
        else
                echo 'LAST_IMAGE' >> $f
        fi

        echo  'ATOMIC_POSITIONS {crystal}' >> $f

##      head -n $d XYZ.dat | tail -n 16  >>  $f

###     Change 1 for S and 2 for H    ###
        head -n $d XYZ.dat | tail -n $n_atoms | awk ' {print "  ", $2, "   ", $3 , "   ", $4 , "   ", $5} ' | sed 's/  1 /  S /g' | sed 's/  2 /  H /g'  >>  $f


        echo $i
done

tail -n 10 $f_in  >> $f

rm XYZ.dat 

echo "##### done #####"


