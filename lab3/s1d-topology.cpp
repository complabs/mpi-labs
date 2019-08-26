/**
    1D simple topology

    # Lab 3 Exercise 2 - Topologies

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
*/

#include <iostream>
#include <iomanip>
#include <mpi.h>

int main( int argc, char** argv )
{
    MPI_Init( &argc, &argv );

    int worldsz = 0, myrank = 0;
    MPI_Comm_size( MPI_COMM_WORLD, &worldsz );
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

    int period = 1;
    MPI_Comm cart_comm;
    MPI_Cart_create( MPI_COMM_WORLD, 1, &worldsz, &period, 0, &cart_comm );

    int cart_id = 0;
    MPI_Comm_rank( cart_comm, &cart_id );

    int cart_position = 0;
    MPI_Cart_coords( cart_comm, myrank, 1, &cart_position );

    int plus_one = 0;
    MPI_Cart_shift( cart_comm, 0, 1, &cart_id, &plus_one );

    int minus_one = 0;
    MPI_Cart_shift( cart_comm, 0, -1, &cart_id, &minus_one );

    std::cout << std::setw(2) << myrank
        << ": cart = " << cart_id
        << ", position = " << cart_position
        << ", plus_one = " << plus_one
        << ", minus_one = " << minus_one << std::endl;

    MPI_Finalize ();

    return 0;
}

