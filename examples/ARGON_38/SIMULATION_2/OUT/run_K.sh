#!/bin/sh
#PJM --rsc-list "node=80"
#PJM --rsc-list "elapse=16:00:00"
#PJM --rsc-list="rscgrp=small"
#PJM --stg-transfiles all
##PJM --mpi "use-rankdir"
#PJM --mpi "proc=640"
#PJM --stgin  "./a.out             ./"
#PJM --stgin  "./atoms.dat    ./"
#PJM --stgin  "./param_3N.in      ./"
##PJM --stgin  "rank=* ./a.out             %r:./"
##PJM --stgin  "rank=* ./atoms.dat    %r:./"
##PJM --stgin  "rank=* ./param_3N.in      %r:./"
##
#PJM --stgout "./*.lammpstrj      ./ "
#PJM --stgout "./coord*.lammpstrj      ./ "
#PJM --stgout "./atomic*.lammpstrj      ./ "
#PJM --stgout "./*.dat     ./ "
#PJM --stgout "./E_atom*.dat     ./ "
#PJM --stgout "./path.dat     ./ "
##PJM --stgout "rank=0 %r:./*.lammpstrj      ./ "
##PJM --stgout "rank=0 %r:./*.dat     ./ "
##PJM --stgout "./out       ./ "
#PJM -s

#---------------------------------------------------------------------------
# path
. /work/system/Env_base_1.2.0-24

# Parameters
NPROCS=8
NTHREADS=1
export PARALLEL=${NPROCS}
export OMP_NUM_THREADS=${NTHREADS}
#LPG="/opt/FJSVxosmmm/sbin/lpgparm -t 4MB -s 4MB -h 4MB -d 4MB -p 4MB"

# Endian
export FORT90L='-Wl,-T'

# MCA Parameters
#MCA_PARAM="--mca common_tofu_fastmode_threshold 0"
#MCA_PARAM="${MCA_PARAM} --mca common_tofu_max_fastmode_procs 40 --mca common_tofu_large_recv_buf_size 2097152"

#---------------------------------------------------------------------------
pwd

date
#./a.out < param_3N.in  > out
mpiexec   ./a.out < param_3N.in  > out
date

