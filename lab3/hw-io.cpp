/**
    Parallel I/O

    # Lab 3 Exercise 3 - MPI I/O

    MPI I/O is used so that results can be written to the same file in
    parallel. Take the serial hello world programs and modify them so that
    instead of writing the output to screen the output is written to a file
    using MPI I/O.

    The simplest solution is likely to be for you to create a character
    buffer, and then use the MPI_File_write_at function.
*/

#include <sstream>
#include <iomanip>
#include <string>
#include <mpi.h>

int main( int argc, char** argv )
{
    MPI_Init( &argc, &argv );

    int myrank, worldsz;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&worldsz);

    // Open the file
    //
    MPI_File fp;
    const char* filename = "result.txt";

    MPI_File_open( MPI_COMM_WORLD,
                   filename, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                   MPI_INFO_NULL, &fp );

    // Write the text into a string stream
    //
    std::ostringstream msg;
    msg << "Rank " << std::setw(3) << myrank << " Hello World!" << std::endl;

    // Copy the string stream result to a buffer
    //
    std::string buf = msg.str();
    MPI_Offset offset = myrank * buf.size();

    MPI_Status status;
    MPI_File_write_at( fp, offset, buf.c_str(), buf.size(), MPI_CHAR, &status );

    // Close the file
    //
    MPI_File_close( &fp );

    MPI_Finalize ();

    return 0;
}
