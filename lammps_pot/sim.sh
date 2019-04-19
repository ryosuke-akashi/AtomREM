#!/bin/sh
#PBS -q gpgpu
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
for node in ${used_nodes}
do
#
    ssh ${node} "
cd ${PBS_O_WORKDIR}
#
# Make local directories in each nodes
#
mkdir -v /scratch/${PBS_JOBID}
#
# (1) Copy-in files to all nodes
#
cp -r -f * /scratch/${PBS_JOBID}
"
#
done
#
# (2) Copy-in files anly to the root node
#
cd ${PBS_O_WORKDIR}
#cp -r xxx.in /scratch/${PBS_JOBID}
#
cd /scratch/${PBS_JOBID}
#
# (3) Program(s) execution
#
mpirun -machinefile ${PBS_NODEFILE} -ppn ${proc_per_node} -np ${num_proc} ${PBS_O_WORKDIR}/a.out < ${PBS_O_WORKDIR}/param_3N.in > out
#
# (4) Copy-out files only from root node
#
cp -r -f * ${PBS_O_WORKDIR}
#
cd ${PBS_O_WORKDIR}
#
for node in ${used_nodes}
do
#
    ssh ${node} "
cd /scratch/${PBS_JOBID} 
#
# (5) Copy-out files from all nodes
#
#cp -r XXX.out ${PBS_O_WORKDIR}
#
# Delete scratch directories
#
cd ${PBS_O_WORKDIR}
rm -r -f /scratch/${PBS_JOBID}
"
#
done
#
echo  
echo      END DATE : `date`
echo  
echo "##### done #####"
