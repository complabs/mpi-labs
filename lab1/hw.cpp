/**
    MPI Hello World

    # Lab 1 Exercise 1: Run "Hello, World"

    Compile and run the "Hello, World" program found in the lecture.
    Make sure you understand how each processors prints its rank as well
    as the total number of processors in the communicator MPI_COMM_WORLD.
*/

#include <iostream>
#include <mpi.h>

int main( int argc, char** argv )
{
    MPI_Init( &argc, &argv );

    int myrank, worldsz;
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    MPI_Comm_size( MPI_COMM_WORLD, &worldsz );

    std::cout << "Processor " << myrank << " of " << worldsz
        << ": Hello World!" << std::endl;

    MPI_Finalize ();

    return 0;
}
