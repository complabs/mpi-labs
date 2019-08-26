/**
    Conway's Game of Life

    # Lab 1 Exercise 5: Use P2P in Game of Life

    In this exercise, you continue learning about point-to-point
    message-passing routines in MPI. After completing this exercise, you
    should be able to write the real parallel MPI code to solve the
    Game of Life.

    To start this exercise, add the initialization and finalization routines
    to the serial "Game of Life" code. This will effectly duplicate the exact
    same calculation on each processor. In order to show that the code is
    performing as expected, add statements to print overall size, and the
    rank of the local process. Don't forget to add the MPI header file.

    ## Domain Decomposition ##

    In order to truly run the "Game of Life" program in parallel, we must
    set up our domain decomposition, i.e., divide the domain into chunks and
    send one chunk to each processor. In the current exercise, we will limit
    ourselves to two processors. If you are writing your code in C, divide
    the domain with a horizontal line, so the upper half will be processed
    on one processor and the lower half on a different processor. If you are
    using Fortran, divide the domain with a vertical line, so the left half
    goes to one processor and the right half to another.

    Hint: Although this can be done with different kinds of sends and receives,
    use blocking sends and receives for the current problem. We have chosen the
    configuration described above because in C arrays, rows are contiguous, and
    in Fortran columns are contiguous. This approach allows the specification
    of the initial array location and the number of words in the send and receive
    routines.

    One issue that you need to consider is that of internal domain boundaries.
    Figure 1 shows the "left-right" domain decomposition described above. Each
    cell needs information from all adjacent cells to determine its new state.
    With domain decomposition, some of the required cells no longer are available
    on the local processor. A common way to tackle this problem is through the use
    of ghost cells. In the current example, a column of ghost cells is added to
    the right side of the left domain, and a column is also added to the left
    side of the right domain (shown in Figure 2). After each time step, the ghost
    cells are filled by passing the appropriate data from the other processor.
    You may want to refer to the figure in the
                [background on the "Game of Life"](Game_of_life.md)
    to see how to fill the other ghost cells.

    Figure 1. Left-right domain decomposition.

    <img src="lr_decomp.jpg" alt="Figure 1" width="400px"/>

    Figure 2. Ghost cells.

    <img src="ghost.jpg" alt="Figure 2" width="400px"/>

    ## Challenge ##

    Implement the domain decomposition described above, and add message passing to
    the ghost cells. Don't forget to divide the domain using a horizontal line for
    C and a vertical line for Fortran. In a subsequent lesson we will examine domain
    decomposition in the opposite direction.
*/

#include "../lab1/golife.h"

///////////////////////////////////////////////////////////////////////////////
/** Implements MPI P2P version of Conway's Game of Life
 */
class GameOfLife_Lab1 : public GameOfLife_MPI
{
public:

    /** Constructs the chunk (allocates memory and initialized arrays).
     * The domain decomposition is in i-direction
     */
    GameOfLife_Lab1
    (
        int mpiProcs, int mpiRank,
        int ni, int nj, bool debugFlag
    )
        : GameOfLife_MPI( mpiProcs, mpiRank, ni / mpiProcs, nj, debugFlag )
    {
    }

protected: //////////////////// VIRTUAL MEMBERS //////////////////////////

    /** Initializes elements of oldc to 0 or 1.
     *  Populate only the local i's on the current rank.
     */
    virtual void InitializeCells ()
    {
        const int offset = myrank * NI;
        for( int i = 1; i <= worldsz * NI; ++i )
        {
            for( int j = 1; j <= NJ; ++j )
            {
                // rand() must be here in the loop to give the same values
                // across all the ranks.
                //
                float x = rand() / ( float(RAND_MAX) + 1);

                if( i >= offset && i <= offset + NI + 1 && x >= 0.5 )
                {
                    oldc[i - offset][j] = 1;
                }
            }
        }
    }

    /** Exchanges the top-bottom boundary of the chunk.
     */
    void exchangeBoundary_TopBottom ()
    {
        if ( worldsz == 1 ) {
            GameOfLife::exchangeBoundary_TopBottom ();
            return;
        }

        int sz = NJ + 2; // Row size
        MPI_Request waitL, waitR;

        //////////////////////////////////////////////////////////////////////////////////
        // Send the bottom row
        //////////////////////////////////////////////////////////////////////////////////
        {
            int ID = myrank * 10; // Tag ID for the message
            int sendTo = myrank >= worldsz - 1 ? 0 : myrank + 1;

            if ( debug ) {
                std::cout << stepNo << " [" << std::setw(2) << myrank << "] " << "Send  ";
                for( int k = 0; k < myrank + 1; k++ ) std::cout << "      ";
                std::cout << std::setw(2) << myrank << " >> " << sendTo << std::endl;
            }

            MPI_Isend( &oldc[NI][0], sz, MPI_INT, sendTo, ID+1, MPI_COMM_WORLD, &waitR );
        }
        //////////////////////////////////////////////////////////////////////////////////
        // Send the top row
        //////////////////////////////////////////////////////////////////////////////////
        {
            int sendTo = myrank <= 0 ? worldsz - 1 : myrank - 1;
            int ID = sendTo* 10; // Tag ID for the message

            if ( debug ) {
                std::cout << stepNo << " [" << std::setw(2) << myrank << "] " << "Send  ";
                for( int k = 0; k < sendTo + 1; k++ ) std::cout << "      ";
                std::cout << std::setw(2) << sendTo << " << " << myrank << std::endl;
            }

            MPI_Isend( &oldc[1][0], sz, MPI_INT, sendTo, ID+2, MPI_COMM_WORLD, &waitL );
        }

        //////////////////////////////////////////////////////////////////////////////////
        // Receive the bottom row (as the ghost)
        //////////////////////////////////////////////////////////////////////////////////
        {
            int ID = myrank * 10; // Tag ID for the message
            int from = myrank >= worldsz - 1 ? 0 : myrank + 1;

            if ( debug ) {
                std::cout << stepNo << " [" << std::setw(2) << myrank << "] " << "Recv  ";
                for( int k = 0; k < myrank + 1; k++ ) std::cout << "      ";
                std::cout << std::setw(2) << myrank << " << " << from << std::endl;
            }

            MPI_Status status;
            MPI_Recv( &oldc[NI+1][0], sz, MPI_INT, from, ID+2, MPI_COMM_WORLD, &status );
        }
        //////////////////////////////////////////////////////////////////////////////////
        // Receive the top row (as the ghost)
        //////////////////////////////////////////////////////////////////////////////////
        {
            int from = myrank <= 0 ? worldsz - 1 : myrank - 1;
            int ID = from * 10; // Tag ID for the message

            if ( debug ) {
                std::cout << stepNo << " [" << std::setw(2) << myrank << "] " << "Recv  ";
                for( int k = 0; k < myrank; k++ ) std::cout << "      ";
                std::cout << std::setw(2) << from << " >> " << myrank << std::endl;
            }

            MPI_Status status;
            MPI_Recv( &oldc[0][0], sz, MPI_INT, from, ID+1, MPI_COMM_WORLD, &status );
        }

        //////////////////////////////////////////////////////////////////////////////////
        // Synchronize with the others
        //
        MPI_Status status;
        MPI_Wait( &waitL, &status );
        MPI_Wait( &waitR, &status );

        if( debug ) { // To get stepNo grouped
            MPI_Barrier( MPI_COMM_WORLD );
        }
    }
};

//////////////////////////////////////////////////////////////////////////////////////////
/** Main program entry.
 *  Initializes MPI, constructs an instance of GameOfLife using
 *  command line parameters, and finally plays the game.
 */
int main( int argc, char** argv )
{
    int total_NI = argc >= 2 ? atoi( argv[1] ) : 200;
    int total_NJ = argc >= 3 ? atoi( argv[2] ) : 200;
    int NSTEPS   = argc >= 4 ? atoi( argv[3] ) : 500;
    bool debug   = argc >= 5 ? atoi( argv[4] ) : false;

    MPI_Init( &argc, &argv );

    int worldsz, myrank;
    MPI_Comm_size( MPI_COMM_WORLD, &worldsz );
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

    if( myrank == 0 )
    {
        std::cout << std::endl
            << "Lab3 Game of Life: NI = " << total_NI << ", NJ = " << total_NJ
            << ", NSTEPS = " << NSTEPS
            << ", MPI world size = " << worldsz << std::endl;
    }

    if( ( total_NI % worldsz ) != 0 )
    {
        if ( myrank == 0 ) {
            std::cerr << "NI must be a multiple of the MPI world size"<< std::endl;
        }
        MPI_Abort( MPI_COMM_WORLD, 1 );
        MPI_Finalize ();
        return -1;
    }

    GameOfLife_Lab1( worldsz, myrank, total_NI, total_NJ, debug )
        .Iterate( NSTEPS );

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize ();

    return 0;
}
