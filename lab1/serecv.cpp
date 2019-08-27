/**
    Send/receive with forwarding

    # Lab 1 Exercise 2: Send data across all processes

    Write a program that takes data from process zero and sends it
    to all of the other processes. That is, process i should receive
    the data and send it to process i+1, until the last process is reached.

    Assume that the data consists of a single integer. For simplicity set
    the value for the first process directly in the code. You may want to
    use MPI_Send and MPI_Recv in your solution.
*/

#include <iostream>
#include <iomanip>
#include <mpi.h>

std::ostream& indent( std::ostream& os, int myrank )
{
    for( int i = 0; i < myrank; ++i ) {
        os << "    ";
    }
    os << std::setw(2) << myrank << ": ";
    return std::cout;
}

int main( int argc, char** argv )
{
    MPI_Init( &argc, &argv );

    int myrank, worldsz;
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    MPI_Comm_size( MPI_COMM_WORLD, &worldsz );

    if( myrank == 0 )
    {
        int value = 77 * 1000;
        MPI_Send( &value, 1, MPI_INT, myrank + 1, 0, MPI_COMM_WORLD );
        indent( std::cout, myrank ) << "Send " << value << std::endl;
    }
    else
    {
        int value = -1;
        MPI_Status status;
        MPI_Recv( &value, 1, MPI_INT, myrank - 1, 0, MPI_COMM_WORLD, &status );
        indent( std::cout, myrank ) << "Recv " << value << std::endl;

        if( myrank < worldsz - 1 ) {
            ++value;
            MPI_Send( &value, 1, MPI_INT, myrank + 1, 0, MPI_COMM_WORLD );
            indent( std::cout, myrank ) << "Send " << value << std::endl;
        }
    }

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize ();

    return 0;
}
