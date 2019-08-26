# Overview

In this lab you will get more familiar with more advanced MPI topics,
including one sided communication and MPI I/O.

# Source Codes

- MPI One sided:    game_of_life-one_sided-prototype.c
- MPI Topology:     simple_1d_topology.c
- MPI I/O:          hello_mpi.c
- MPI Latency:      mpi_latency.c
- MPI Bandwidth:    mpi_bandwidth.c
- MPI Bandwidth NB: mpi_bandwidth-nonblock.c

----------------------------------------------------------------------------------

# Exercise 1 - One sided communication

Take the prototype one sided communication code and complete the code by
adding the correct one sided MPI calls so that the program works. The
number of live cells after the calculation should be the same on any
number of tasks that can easily divide the grid. The solution that will
be provided towards the end of the class uses MPI_Get, but something
similar could also be done with MPI_Put.

----------------------------------------------------------------------------------

# Exercise 2 - Topologies

### Part A Run the simple example topology program and understand how it
works. Notice that the rank order in the MPI_COMM_WORLD communicator is
not necessarily the same as for the cart_comm communicator.

### Part B The code in Exercise 1 uses a simple and manually implemented
"topology". Re-implement the calculation of which MPI task to read the
halo data from using MPI topology functions, i.e. set up a simple
periodic 1d topology then use MPI_Cart_shift to get the rank of the
ranks to get the data from.

Note that the position in the new topology is not necessarily the same
as the position in MPI_COMM_WORLD so make sure that the initial grid
setup reflects that.

----------------------------------------------------------------------------------

# Exercise 3 - MPI I/O

MPI I/O is used so that results can be written to the same file in
parallel. Take the serial hello world programs and modify them so that
instead of writing the output to screen the output is written to a file
using MPI I/O.

The simplest solution is likely to be for you to create a character
buffer, and then use the MPI_File_write_at function.

----------------------------------------------------------------------------------

# Exercises 4 - Bandwidth and latency between nodes

Use `mpi_wtime` to compute latency and bandwidth with the bandwidth and
latency codes above

**Note**: In modifying the original exercises provided by LLNL, We had
to make a small change to the latency code as the Cray latency is a lot
better than the tests were designed for. When the latency is of the
order 1 millisecond, writing it out as an integer number of milliseconds
did not make much sense.

For this exercise, it is nice to compare running on the same node e.g.

``` salloc -N 1 --ntasks-per-node=2 -A <project> -t 00:05:00 srun -n 2
./mpi_latency.x ```

with running on separate nodes

``` salloc -N 2 --ntasks-per-node=1 -A <project> -t 00:05:00 srun -n 2
./mpi_latency.x ```

Similarly for the bandwidth.

As you would expect the latency is much better on a single node than
across nodes, but possibly unexpectedly if you just have 2 MPI tasks the
bandwidth is better between nodes than across a single node. (probably
related to lack of contention for resources, e.g. the gemini chips and
the l3 cache etc.)

----------------------------------------------------------------------------------

# Acknowledgment

The examples in this lab are provided for educational purposes by
[National Center for Supercomputing
Applications](http://www.ncsa.illinois.edu/), (in particular their
[Cyberinfrastructure Tutor](http://www.citutor.org/)), [Lawrence
Livermore National Laboratory](https://computing.llnl.gov/) and [Argonne
National Laboratory](http://www.mcs.anl.gov/). Much of the LLNL MPI
materials comes from the [Cornell Theory
Center](http://www.cac.cornell.edu/). We would like to thank them for
allowing us to develop the material for machines at PDC. You might find
other useful educational materials at these sites.

