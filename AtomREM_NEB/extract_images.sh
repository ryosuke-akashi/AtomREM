#!/bin/sh
#

f='inout_AtomREM/atomic_coord_q.lammpstrj' 
echo $f

f2='inout_AtomREM/coordinates_rank_100020.lammpstrj'
echo $f2

num=9

##last_step=`grep -A1 'ITEM: TIMESTEP'  $f | tail -n 1`

##pr_step=`grep -A1 'ITEM: TIMESTEP'  $f | grep "[0-9]" | tail -n 2 | head -n 1 `

##time_step=`expr $last_step - $pr_step`
##echo $time_step
 

##first_step=`grep -A1 'ITEM: TIMESTEP'  $f | head -n 2 | tail -n 1 `


##num_steps=`grep -A 1 'ITEM: TIMESTEP'  $f | grep -v 'ITEM: TIMESTEP' | grep "[0-9]" | wc -w`
#echo $num_steps

n_atoms=`head -n 10 $f | grep -A 1 'ITEM: NUMBER OF ATOMS'  | grep -v 'ITEM: NUMBER OF ATOMS' | grep "[0-9]"`
n_lines=`expr $n_atoms + 9`

num_im=`wc -l < $f`
num_im=`expr $num_im \/ $n_lines`
echo $num_im

step=`expr $num_im / $num` 
echo $step

rm images.dat 

for ((i=1; i <= $num_im - $step ; i = i + $step ))
do


echo $i

st=`expr  $i  \* $n_lines` 
echo $st

head -n ${st}   $f  | tail -n $n_lines  >>  images.dat

done


nall_1=`expr ${num} \* ${n_lines} `

echo $nall_1

head -n ${nall_1}   images.dat > im_2.dat 

tail -n $n_lines $f  >> im_2.dat 



mv im_2.dat images.dat 



###################  ADD the configurations from relaxation path  ################################

f=$f2
echo $f


num=7

num_im=`wc -l < $f`
num_im=`expr $num_im \/ $n_lines`
echo $num_im

step=`expr $num_im / $num`
echo $step

for ((i=1; i <= $num_im - $step ; i = i + $step ))
do


echo $i

st=`expr  $i  \* $n_lines`
echo $st

head -n ${st}   $f  | tail -n $n_lines  >>  images_2.dat

done


nall_1=`expr ${num} \* ${n_lines} `

echo $nall_1

head -n ${nall_1}   images_2.dat > im_2.dat

tail -n $n_lines $f  >> im_2.dat



#################  MERGE


mv im_2.dat images_2.dat 

cat images_2.dat >> images.dat	


rm images_2.dat


