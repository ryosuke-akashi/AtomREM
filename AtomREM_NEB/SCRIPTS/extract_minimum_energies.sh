#!/bin/sh
#
 echo " extract minimal energy after PW relaxation for each image:   sh extract_minimum_energies.sh "


date



grep -n "! "  IM_*/H3S.out  > 1.txt

sed  -e 's/IM_/ /g' 1.txt > 2.txt

sed  -e 's/\/H3S/     /g' 2.txt > 3.txt

awk ' {print "  ", $1, "   ", $8} ' 3.txt | sort -n   >  path_relax.dat

sort -n -r path_relax.dat | sort -n -u -k1  > 6.txt

paste  H3S_NEB.dat   6.txt  | awk ' {print $4, "   ", $1, "     ", $5} '   >    path_min.dat


rm ex.txt

while IFS= read -r line; do

    num=`echo $line | cut -d " " -f  1 `
    rc=`echo $line | cut -d " " -f  2 `

    echo $num "   "  $rc

    grep -n  " ${num}  "    path_relax.dat     |    sed  -e 's/'"${num}  "'/'"${rc}"'/g'   >> ex.txt

done   <     path_min.dat


awk ' {print "  ", $2, "   ", $3} ' ex.txt    >  path_rc_relax.dat


rm {1,2,3,6}.txt
rm ex.txt

echo "##### done #####"


