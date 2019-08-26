#!/bin/bash -l

#SBATCH -J mpi-L1

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
    SRUN=mpirun
fi


echo "/////////////////////////////////////////////////////////////////////////"
echo "//                              LAB  1                                 //"
echo "/////////////////////////////////////////////////////////////////////////"

echo
echo "============ Lab 1 Exercise 1: Run Hello, World"
echo

$SRUN -n 32 ./hw

echo
echo "============ Lab 1 Exercise 2: Send data across all processes"

echo

$SRUN -n 32 ./serecv

echo
echo "============ Lab 1 Exercise 3: Find pi using P2P Communication (Master/Worker)"
echo

$SRUN -n 32 ./pi

echo
echo "============ Now run with 10 workers"
echo

$SRUN -n 10 ./pi

echo
echo "============ Lab 1 Exercise 4: Use P2P in Parallel Search"
echo

$SRUN -n 32 ./psearch

echo
echo "============ Lab 1 Exercise 5: Use P2P in Game of Life"
echo

$SRUN -n 1 ./golife # Reference
$SRUN -n 2 ./golife
$SRUN -n 10 ./golife

echo

