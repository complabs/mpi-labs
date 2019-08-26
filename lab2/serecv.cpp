/**
    Send/receive with forwarding

    # Lab 2 Exercise 2: Send data across all processes using non-blocking

    Take the code for sending data across all processes from the MPI Lab 1,
    and have each node add one to the number received, print out the result,
    and send the results on.

    ### Use Proper Synchronization

    For the case where you want to use proper synchronization, you'll want
    to do a non-blocking receive, add one, print, then a non-blocking send.
    The result should be `1 - 2 - 3 - 4 - 5 ...`

    ### Try without Synchronization: Detect Race Conditions

    To see what happens without synchronization, leave out the `wait`.
*/

#include <iostream>
#include <mpi.h>

int main( int argc, char** argv )
{
    /// A race condition flag (if given any argument introduce a race condition)
    bool race = argc >= 2;

    MPI_Init( &argc, &argv );

    int myrank, worldsz;
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    MPI_Comm_size( MPI_COMM_WORLD, &worldsz );

    if( myrank == 0 ) // master = the first
    {
        int value = 77 * 1000;
        std::cout << myrank << ": Sending " << value << std::endl;
        MPI_Send( &value, 1, MPI_INT, myrank + 1, 0, MPI_COMM_WORLD );

        MPI_Request request;
        MPI_Isend( &value, 1, MPI_INT, myrank + 1, 0, MPI_COMM_WORLD,&request );
    }
    else // ranks in the chain
    {
        int value = -1;
        MPI_Request request;
        MPI_Irecv( &value, 1, MPI_INT, myrank - 1, 0, MPI_COMM_WORLD, &request );

        if ( race )
        {
            std::cout << myrank << ": After MPI_Irecv (nowait) " << value << std::endl;
        }
        else
        {
            MPI_Status status;
            MPI_Wait( &request, &status );

            std::cout << myrank << ": Received " << value << std::endl;
        }

        ++value;
        if( myrank < worldsz - 1 )
        {
            std::cout << myrank << ": Forwarded " << value << std::endl;
            MPI_Request request2;
            MPI_Isend( &value, 1, MPI_INT, myrank + 1, 0, MPI_COMM_WORLD, &request2 );
        }
    }

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize ();

    return 0;
}
