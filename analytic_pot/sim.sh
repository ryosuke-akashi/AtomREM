#!/bin/sh
#PBS -q has2667 
#PBS -j oe
#PBS -l nodes=1:ppn=16
#
# You should edit following sections before you submit job.
#
# (1) Copy-in files to all nodes
# (2) Copy-in files only to root node
# (3) Program(s) execution
# (4) Copy-out files only from root node
# (5) Copy-out files from all nodes
#
#############  Initialization (You should not modify here)  #############
# 
source /opt/intel/parallel_studio_xe_2015/psxevars.sh
#
proc_per_node=1
#
export OMP_NUM_THREADS=`expr 16 / ${proc_per_node}`
used_nodes=`awk 'NR%16==1{print $1}' ${PBS_NODEFILE}`
num_nodes=`wc -l ${PBS_NODEFILE} | awk '{print $1}'`
#num_nodes=`expr ${num_nodes} / 16`
num_proc=`expr ${proc_per_node} \* ${num_nodes}`
#
echo 
echo    START DATE : `date`
echo     PBS_JOBID : ${PBS_JOBID}
echo PBS_O_WORKDIR : ${PBS_O_WORKDIR}
echo   PBS_NODENUM : ${PBS_NODENUM}
echo    PBS_O_HOST : ${PBS_O_HOST}
echo    PBS_SERVER : ${PBS_SERVER}
echo   PBS_O_QUEUE : ${PBS_O_QUEUE}
echo   PBS_ARRAYID : ${PBS_ARRAYID}
echo     PBS_QUEUE : ${PBS_QUEUE}
echo      used_nodes : ${used_nodes}
echo       num_nodes : ${num_nodes}
echo   proc_per_node : ${proc_per_node}
echo        num_proc : ${num_proc}
echo OMP_NUM_THREADS : ${OMP_NUM_THREADS}
echo     "hostname : " `hostname`
#
#######################  End of Initialization   ##############################
#                                                                                
#
# (3) Program(s) execution
#
cd ${PBS_O_WORKDIR}
LD_LIBRARY_PATH="/home/iurii/src/lammps-5Jun19/src/:"$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

mpirun -machinefile ${PBS_NODEFILE} -ppn ${proc_per_node} -np ${num_proc} ./a.out < ./param_3N.in > ./out
#
echo      END DATE : `date`
echo  
echo "##### done #####"
