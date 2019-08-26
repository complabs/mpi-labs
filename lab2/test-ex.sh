#!/bin/bash -l

#SBATCH -J mpi-L2

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
echo "//                              LAB  2                                 //"
echo "/////////////////////////////////////////////////////////////////////////"

echo
echo "============ Lab 2 Exercise 1: Calculating pi using collectives"
echo

$SRUN -n 32 ./pi

echo
echo "============ Lab 2 Exercise 2.1: Sending data using non-blocking"
echo

$SRUN -n 32 ./serecv

echo
echo "============ Lab 2 Exercise 2.2: Sending data using non-blocking with race"
echo

$SRUN -n 32 ./serecv race

echo
echo "============ Lab 2 Exercise 3: Finding pi using non-blocking communications"
echo

$SRUN -n 32 ./pi-nb

echo
echo "============ Lab 2 Exercise 4.1: Parallel search using collectives"
echo

$SRUN -n 30 ./psearch

echo
echo "============ Lab 2 Exercise 4.1: Game of Life using collectives"
echo

$SRUN -n 1 ./golife  # Reference
$SRUN -n 2 ./golife
$SRUN -n 10 ./golife

# for x in `seq 0 1 3`; do
    # echo "File life-$x.gol -----------------------------------------------------------------------------"
    # cat life-$x.gol
# done

echo
