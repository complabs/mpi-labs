/**
    Parallel search

    # Lab 1 Exercise 4: Use P2P in "Parallel Search"

    In this exercise, you learn about the heart of MPI:
    point-to-point message-passing routines in both their blocking
    and non-blocking forms as well as the various modes of communication.

    Your task is to parallelize the "Parallel Search" problem.
    In the parallel search problem, the program should find all occurrences
    of a certain integer, which will be called the target. It should then
    write the target value, the indices and the number of occurences to
    an output file. In addition, the program should read both the target
    value and all the array elements from an input file.

    Hint: One issue that comes up when parallelizing a serial code is
    handling I/O. As you can imagine, having multiple processes writing to
    the same file at the same time can produce useless results. A simple
    solution is to have each process write to an output file named with its
    rank. Output to these separate files removes the problem.
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

    b_data( const std::string filename, int size )
        : target( 0 ), len( size ), b( NULL )
    {
        std::ifstream infile( filename.c_str () );

        infile >> target;

        b = new int[ len ];
        for( int i = 0; i < len; ++i )
        {
            infile >> b[i];
        }
    }

    ~b_data ()
    {
        delete[] b;
    }
};

int main( int argc, char** argv )
{
    MPI_Init( &argc, &argv );

    int myrank, worldsz;
    MPI_Comm_size( MPI_COMM_WORLD, &worldsz );
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

    // A tag received from slaves at the end of work
    //
    const int MSG_QUIT_TAG = 77;
    const int MSG_HEAD_TAG = 78;
    const int MSG_DATA_TAG = 79;
    const int MSG_INDX_TAG = 80;

    if( myrank == 0 ) // Master
    {
        b_data data( "b.data", 300 );

        // Divide b_data into slaves and chunks. The last chunk may be larger
        // than the others. Example:
        //      nr_slaves = 31
        //      chunk_size = 300 / 31 = 9
        //      lastchunk_size = 300 - 30 * 9 = 30
        //
        int nr_slaves = worldsz - 1;
        int chunk_size = data.len / nr_slaves;
        int lastchunk_size = data.len - ( nr_slaves - 1 ) * chunk_size;

        for( int i = 0; i < nr_slaves; ++i )
        {
            // send header: target, offset, length
            //
            int header[ 3 ] = {
                data.target,
                i * chunk_size,
                i == nr_slaves - 1 ? lastchunk_size : chunk_size
            };
            MPI_Send( header, 3, MPI_INT, i+1, MSG_HEAD_TAG, MPI_COMM_WORLD );

            // send data
            //
            MPI_Send( data.b + header[1], header[2],
                      MPI_INT, i+1, MSG_DATA_TAG, MPI_COMM_WORLD);
        }

        int nr_completed = 0;

        while( nr_completed != nr_slaves )
        {
            MPI_Status status;
            int index;
            MPI_Recv( &index, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG,
                     MPI_COMM_WORLD, &status);

            if( status.MPI_TAG == MSG_QUIT_TAG ) {
                ++nr_completed;
            }
            else {
                std::cout << index << std::endl;
            }
        }
    }
    else // Slave
    {
        // Receive the header with target, offset and data length
        //
        int header[3];
        MPI_Status status;
        MPI_Recv( header, 3, MPI_INT, 0, MSG_HEAD_TAG, MPI_COMM_WORLD, &status );

        int  target   = header[0];
        int  b_offset = header[1];
        int  b_len    = header[2];
        int* b_data   = new int[ b_len ];

        // Receive the data
        //
        MPI_Recv( b_data, b_len, MPI_INT, 0, MSG_DATA_TAG, MPI_COMM_WORLD, &status );

        // Process the data
        //
        for( int i = 0; i < b_len; ++i )
        {
            if( b_data[i] == target )
            {
                int index = b_offset + i + 1;
                MPI_Send( &index, 1, MPI_INT, 0, MSG_INDX_TAG, MPI_COMM_WORLD );
            }
        }

        // Indicate that we are quitting
        //
        MPI_Send( &target, 1, MPI_INT, 0, MSG_QUIT_TAG, MPI_COMM_WORLD );

        delete[] b_data;
    }

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize ();

    return 0;
}
