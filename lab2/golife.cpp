/**
    Conway's Game of Life

    # Lab 2 Exercise 4.2: Implement the Game of Life using collectives

    In almost every MPI program there are instances where all the processors
    in a communicator need to perform some sort of data transfer or
    calculation. These "collective communication" routines are the subject
    of this exercise and the "Parallel Search" and "Game of Life" programs
    are no exception.

    Modify your previous "Game of Life" code to use `mpi_reduce` to compute
    the total number of live cells, rather than individual sends and
    receives.
*/

#include "../lab1/golife.h"

///////////////////////////////////////////////////////////////////////////////
/** Implements MPI version of Conway's Game of Life using collectives
 */
class GameOfLife_Lab2 : public GameOfLife_MPI
{
    MPI_Datatype COLUMN_VEC; //!< Derived type for a column vector.

public:

    /** Constructs the chunk (allocates memory and initialized arrays).
     * The domain decomposition is in i-direction
     */
    GameOfLife_Lab2
    (
        int mpiProcs, int mpiRank,
        int ni, int nj, bool debugFlag
    )
        : GameOfLife_MPI( mpiProcs, mpiRank, ni, nj / mpiProcs, debugFlag )
    {
        // We have to exchange the column boundaries,
        // hence we need a column type vector of length NI+2.
        //
        MPI_Type_vector
        (
            NI + 2,  // Number of blocks
            1,       // Number of elements in each block
            NJ + 2,  // Stride  between the blocks
            MPI_INT, // Element type
            &COLUMN_VEC //
        );
        MPI_Type_commit( &COLUMN_VEC );
    }

public: //////////////////// VIRTUAL MEMBERS //////////////////////////

    /** Initializes elements of oldc to 0 or 1.
     *  Populate only the local j's on the current rank.
     */
    virtual GameOfLife& InitializeCells ()
    {
        const int offset = myrank * NJ;
        for( int i = 1; i <= NI; ++i )
        {
            for( int j = 1; j <= worldsz * NJ; ++j )
            {
                // rand() must be here in the loop to give the same values
                // across all the ranks.
                //
                float x = rand() / ( float(RAND_MAX) + 1);

                if( j >= offset && j <= offset + NJ + 1 && x >= 0.5 )
                {
                    oldc[i][j - offset] = 1;
                }
            }
        }

        return *this;
    }

protected: //////////////////// VIRTUAL MEMBERS //////////////////////////

    /** Exchanges the left-right boundary of the chunk.
     */
    virtual void exchangeBoundary_LeftRight ()
    {
        if ( worldsz == 1 ) {
            GameOfLife::exchangeBoundary_LeftRight ();
            return;
        }

        MPI_Request waitL, waitR;

        //////////////////////////////////////////////////////////////////////////////////
        // Send the right column
        //////////////////////////////////////////////////////////////////////////////////
        {
            int ID = myrank * 10; // Tag ID for the message
            int sendTo = myrank >= worldsz - 1 ? 0 : myrank + 1;

            if ( debug ) {
                std::cout << stepNo << " [" << std::setw(2) << myrank << "] " << "Send  ";
                for( int k = 0; k < myrank + 1; k++ ) std::cout << "      ";
                std::cout << std::setw(2) << myrank << " >> " << sendTo << std::endl;
            }

            MPI_Isend( &oldc[0][NJ], 1, COLUMN_VEC, sendTo, ID+1, MPI_COMM_WORLD, &waitR );
        }
        //////////////////////////////////////////////////////////////////////////////////
        // Send the left column
        //////////////////////////////////////////////////////////////////////////////////
        {
            int sendTo = myrank <= 0 ? worldsz - 1 : myrank - 1;
            int ID = sendTo* 10; // Tag ID for the message

            if ( debug ) {
                std::cout << stepNo << " [" << std::setw(2) << myrank << "] " << "Send  ";
                for( int k = 0; k < sendTo + 1; k++ ) std::cout << "      ";
                std::cout << std::setw(2) << sendTo << " << " << myrank << std::endl;
            }

            MPI_Isend( &oldc[0][1], 1, COLUMN_VEC, sendTo, ID+2, MPI_COMM_WORLD, &waitL );
        }

        //////////////////////////////////////////////////////////////////////////////////
        // Receive the right column (as the ghost)
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
            MPI_Recv( &oldc[0][NJ+1], 1, COLUMN_VEC, from, ID+2, MPI_COMM_WORLD, &status );
        }
        //////////////////////////////////////////////////////////////////////////////////
        // Receive the left column (as the ghost)
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
            MPI_Recv( &oldc[0][0], 1, COLUMN_VEC, from, ID+1, MPI_COMM_WORLD, &status );
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
            << "Lab2 Game of Life: NI = " << total_NI << ", NJ = " << total_NJ
            << ", NSTEPS = " << NSTEPS
            << ", MPI world size = " << worldsz << std::endl;
    }

    if( ( total_NJ % worldsz ) != 0 )
    {
        if ( myrank == 0 ) {
            std::cerr << "NI must be a multiple of the MPI world size"<< std::endl;
        }
        MPI_Abort( MPI_COMM_WORLD, 1 );
        MPI_Finalize ();
        return -1;
    }

    GameOfLife_Lab2( worldsz, myrank, total_NI, total_NJ, debug )
        .InitializeCells ()
        .Iterate( NSTEPS )
        .FinalizeEvolution ();

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize ();

    return 0;
}
