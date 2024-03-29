#!/bin/bash -l

#SBATCH -J mpi-L3

#SBATCH -A edu19.summer
# BATCH --reservation=summer-2019-08-23

#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH -t 00:10:00

#SBATCH -o test-ex.out

if [ ".$SLURM_CLUSTER_NAME" = ".tegner" ]; then
    module load gcc openmpi
    SRUN="mpirun --mca orte_base_help_aggregate 0 --mca  btl_openib_allow_ib true"
elif [ ".$SLURM_CLUSTER_NAME" = ".beskow" ]; then
    SRUN=srun
else
    SLURM_CLUSTER_NAME=`hostname -s`
    SRUN=mpirun
fi


echo "/////////////////////////////////////////////////////////////////////////"
echo "//                              LAB  3                                 //"
echo "/////////////////////////////////////////////////////////////////////////"

echo
echo "============ Lab 3 Exercise 1 - One sided communication"
echo

$SRUN -n 10 ./golife-win

echo
echo "============ Lab 3 Exercise 2.1 - Topologies, game of life"
echo

$SRUN -n 10 ./golife-topology

echo
echo "============ Lab 3 Exercise 2.2 - Topologies, simple 1D"
echo

$SRUN -n 32 ./s1d-topology

echo
echo "============ Lab 3 Exercise 3 - MPI I/O (simple)"
echo

$SRUN -n 32 ./hw-io

echo
echo "============ Lab 3 Exercise 4.1 - Bandwidth and latency between nodes"
echo

$SRUN -n 32 bw > test-${SLURM_CLUSTER_NAME}-bw.out

echo
echo "============ Lab 3 Exercise 4.2 - Non-blocking bandwidth"
echo

$SRUN -n 32 bw-nb > test-${SLURM_CLUSTER_NAME}-bwnb.out

echo
echo "============ Lab 3 Exercise 4.3 - Latency"
echo

$SRUN -n 32 latency > test-${SLURM_CLUSTER_NAME}-lat.out

echo
echo "============ Lab 3 New Exercise 1 - MPI I/O (STL file)"
echo

$SRUN -n 1 ./dtypes data/sphere.stl out1.stl
echo
md5sum out1.stl data/sphere.stl

echo 

$SRUN -n 7 ./dtypes data/sphere.stl out2.stl
echo
md5sum out2.stl data/sphere.stl

echo

$SRUN -n 2 ./dtypes data/sphere.stl out2.stl 2 10
echo
$SRUN -n 10 ./dtypes data/sphere.stl out2.stl 2 10
echo
$SRUN -n 20 ./dtypes data/sphere.stl out2.stl 2 10
echo
$SRUN -n 30 ./dtypes data/sphere.stl out2.stl 2 10

echo

echo
