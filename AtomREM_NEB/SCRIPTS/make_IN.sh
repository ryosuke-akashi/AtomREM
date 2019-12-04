#!/bin/sh
#

 echo " run the NEB calcultion with short cutoff:  sh  make_IN.sh  "

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
                head -n 48 $f_in > $f
             elif [ $i -lt 18 ] ; then
                echo 'INTERMEDIATE_IMAGE' >> $f
             else
                echo 'LAST_IMAGE' >> $f
        fi

        echo  'ATOMIC_POSITIONS {crystal}' >> $f

##      head -n $d XYZ.dat | tail -n 16  >>  $f

###     Change 1 for S and 2 for H    ###

        if [ $i -eq 1 ]; then
             grep -A17 'FIRST_IMAGE' $f_in | tail -n 16  >>  $f  
           else 
             head -n $d XYZ.dat | tail -n $n_atoms | awk ' {print "  ", $2, "   ", $3 , "   ", $4 , "   ", $5} ' | sed 's/  1 /  S /g' | sed 's/  2 /  H /g'  >>  $f
        fi

        echo $i
done


echo 'END_POSITIONS' >> $f

grep -A1 'K_POINTS' $f_in >> $f

grep -A3 'CELL_PARAMETERS' $f_in >> $f

echo 'END_ENGINE_INPUT' >> $f
echo 'END' >> $f







### tail -n 10 $f_in  >> $f

rm XYZ.dat 

echo "##### done #####"


