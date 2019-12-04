#!/bin/sh
#

 echo "simulation of PW relaxation: sh run.sh"

date

sh ../../../environment_variables

### mpijob -np 432  pw.x -npool 18   <  H3S.in   >  H3S.out 

 echo " BULK parallel pw relaxation for each image independently: " 

mpijob.bulk  -np 24   pw.x  <  IM_1/H3S.in   >  IM_1/H3S.out    & 
sleep   5 
mpijob.bulk  -np 24   pw.x  <  IM_2/H3S.in   >  IM_2/H3S.out    &
sleep   5
mpijob.bulk  -np 24   pw.x  <  IM_3/H3S.in   >  IM_3/H3S.out    &
sleep   5
mpijob.bulk  -np 24   pw.x  <  IM_4/H3S.in   >  IM_4/H3S.out    &
sleep   5
mpijob.bulk  -np 24   pw.x  <  IM_5/H3S.in   >  IM_5/H3S.out    &
sleep   5
mpijob.bulk  -np 24   pw.x  <  IM_6/H3S.in   >  IM_6/H3S.out    &
sleep   5
mpijob.bulk  -np 24   pw.x  <  IM_7/H3S.in   >  IM_7/H3S.out    &
sleep   5
mpijob.bulk  -np 24   pw.x  <  IM_8/H3S.in   >  IM_8/H3S.out    &
sleep   5
mpijob.bulk  -np 24   pw.x  <  IM_9/H3S.in   >  IM_9/H3S.out    &
sleep   5
mpijob.bulk  -np 24   pw.x  <  IM_10/H3S.in   >  IM_10/H3S.out    &
sleep   5
mpijob.bulk  -np 24   pw.x  <  IM_11/H3S.in   >  IM_11/H3S.out    &
sleep   5
mpijob.bulk  -np 24   pw.x  <  IM_12/H3S.in   >  IM_12/H3S.out    &
sleep   5
mpijob.bulk  -np 24   pw.x  <  IM_13/H3S.in   >  IM_13/H3S.out    &
sleep   5
mpijob.bulk  -np 24   pw.x  <  IM_14/H3S.in   >  IM_14/H3S.out    &
sleep   5
mpijob.bulk  -np 24   pw.x  <  IM_15/H3S.in   >  IM_15/H3S.out    &
sleep   5
mpijob.bulk  -np 24   pw.x  <  IM_16/H3S.in   >  IM_16/H3S.out    &
sleep   5
mpijob.bulk  -np 24   pw.x  <  IM_17/H3S.in   >  IM_17/H3S.out    &
sleep   5
mpijob.bulk  -np 24   pw.x  <  IM_18/H3S.in   >  IM_18/H3S.out    &
wait




echo  
echo "##### done: run.sh  #####"

date

