/**
    Parallel search

    # Lab 2 Exercise 4.1: Implement the "Parallel Search" Using Collectives

    In almost every MPI program there are instances where all the processors
    in a communicator need to perform some sort of data transfer or
    calculation. These "collective communication" routines are the subject
    of this exercise and the "Parallel Search" and "Game of Life" programs
    are no exception.

    ### Task:

    Modify your previous "Parallel Search" code to change how the master
    first sends out the target and subarray data to the slaves. Use the MPI
    broadcast routines to give each slave the target. Use the MPI scatter
    routine to give all processors a section of the array ``b`` it will
    search.

    Hint: When you use the standard MPI scatter routine you will see that
    the global array ``b`` is now split up into four parts and the master
    process now has the first fourth of the array to search. So you should
    add a search loop (similar to the workers') in the master section of
    code to search for the target and calculate the average and then write
    the result to the output file. This is actually an improvement in
    performance since all the processors perform part of the search in
    parallel.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <mpi.h>

// Encapsulates `b.data` file which has the target value on the first line
// and the remaining `len`-lines with the values for the `b` array.
//
class b_data
{
public:
    int target;
    int len;
    int* b;

    b_data( int size )
        : target( 0 ), len( size ), b( NULL )
    {
    }

    void allocate ()
    {
        if ( b == NULL ) {
            b = new int[ len ];
        }
    }

    void load( const std::string filename )
    {
        allocate ();

        std::ifstream infile( filename.c_str () );

        infile >> target;

        for( int i = 0; i < len; ++i ) {
            infile >> b[i];
        }
    }

    ~b_data ()
    {
        if ( b != NULL ) {
            delete[] b;
        }
    }
};

int main( int argc, char** argv )
{
    MPI_Init( &argc, &argv );

    int myrank, worldsz;
    MPI_Comm_size( MPI_COMM_WORLD, &worldsz );
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

    // Setup the problem size
    //
    b_data data( 300 );

    if( ( data.len % worldsz ) != 0 )
    {
        if ( myrank == 0 ) {
            std::cerr << "N must be a multiple of the MPI world size"<< std::endl;
        }
        MPI_Abort( MPI_COMM_WORLD, 1 );
        MPI_Finalize ();
        return -1;
    }

    // Master can now load the data
    //
    if( myrank == 0 ) {
        data.load( "b.data" );
    }

    // Send the target (called by all ranks)
    //
    MPI_Bcast( &data.target, 1, MPI_INTEGER, 0, MPI_COMM_WORLD );

    // The data is scattered, here is our local chunk
    //
    b_data local( data.len / worldsz );
    local.allocate ();

    MPI_Scatter( data.b, local.len, MPI_INTEGER, local.b, local.len, MPI_INTEGER,
                 0, MPI_COMM_WORLD );

    // Search the b array and save the target locations and number
    //
    b_data result( local.len );
    result.allocate ();
    int result_count = 0;

    for( int i = 0; i < local.len; ++i )
    {
        // If we found the target, save its offset as a result
        //
        if( local.b[i] == data.target ) {
            result.b[ result_count++ ] = myrank * local.len + i + 1;
        }
    }

    // Gather the partial count from each process
    //
    int partial_count[ worldsz ];
    MPI_Gather( &result_count, 1, MPI_INTEGER, partial_count, 1, MPI_INTEGER,
                0, MPI_COMM_WORLD );

    // Note: to avoid waisting the memory, the results are sored into data.b
    //
    if( myrank == 0 ) // Master calculates the total count and displays the results
    {
        int global_len = 0;
        int displs[ worldsz ]; // the displacement relative to global_res
        for( int i = 0; i < worldsz; ++i )
        {
            displs[i] = global_len;
            global_len += partial_count[i];
        }

        MPI_Gatherv( result.b, result_count, MPI_INTEGER,
                     /* recvbuf */ data.b, /* recvcounts */ partial_count,
                     displs, MPI_INTEGER, 0, MPI_COMM_WORLD );

        for( int i = 0; i < global_len; ++i ) {
            std::cout << data.b[i] << std::endl;
        }
    }
    else // Slave just gathers into specified locations
    {
        int displs[ worldsz ]; // the displacement relative to global_res
        MPI_Gatherv( result.b, result_count, MPI_INTEGER,
                     /* recvbuf */ data.b, /* recvcounts */ partial_count,
                     displs, MPI_INTEGER, 0, MPI_COMM_WORLD );
    }

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize ();

    return 0;
}
