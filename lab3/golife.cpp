/**
    Conway's Game of Life

    # Lab 3 Exercise 1 - One sided communication

    Take the prototype one sided communication code and complete the code by
    adding the correct one sided MPI calls so that the program works. The
    number of live cells after the calculation should be the same on any
    number of tasks that can easily divide the grid. The solution that will
    be provided towards the end of the class uses MPI_Get, but something
    similar could also be done with MPI_Put.

    # Lab 3 Exercise 2 - Topologies

    To enable: #define TOPOLOGY

    ## Part A

    Run the simple example topology program and understand how it
    works. Notice that the rank order in the MPI_COMM_WORLD communicator is
    not necessarily the same as for the cart_comm communicator.

    ## Part B

    The code in Exercise 1 uses a simple and manually implemented
    "topology". Re-implement the calculation of which MPI task to read the
    halo data from using MPI topology functions, i.e. set up a simple
    periodic 1d topology then use MPI_Cart_shift to get the rank of the
    ranks to get the data from.

    Note that the position in the new topology is not necessarily the same
    as the position in MPI_COMM_WORLD so make sure that the initial grid
    setup reflects that.
*/

#include "../lab1/golife.h"

///////////////////////////////////////////////////////////////////////////////
/** Implements Conway's Game of Life with one sided MPI communication
 */
class GameOfLife_Lab3 : public GameOfLife_MPI
{
    MPI_Win top_win;     //!< One-sided window to the top row cells
    MPI_Win bottom_win;  //!< One-sided window to the bottom row cells
    int my_position;

#ifdef TOPOLOGY
    MPI_Comm cart_comm;  //!< An instance of the periodic 1D topology
    int cart_id;         //!< Our identifier in the topology
#endif

public:

    /** Constructs the chunk (allocates memory and initialized arrays).
     * The domain decomposition is in i-direction
     */
    GameOfLife_Lab3
    (
        int mpiProcs, int mpiRank,
        int ni, int nj, bool debugFlag
    )
        : GameOfLife_MPI( mpiProcs, mpiRank, ni / mpiProcs, nj, debugFlag )
        , my_position( mpiRank )
    {
        // Create the memory windows to our top/bottom rows
        //
        int dunits = sizeof(int);
        MPI_Aint sz = ( NJ + 2 ) * sizeof(int);
        MPI_Win_create( &oldc[1][0],  sz, dunits, MPI_INFO_NULL, MPI_COMM_WORLD, &top_win    );
        MPI_Win_create( &oldc[NI][0], sz, dunits, MPI_INFO_NULL, MPI_COMM_WORLD, &bottom_win );

#ifdef TOPOLOGY
        // Create a periodic 1D topology.
        // Set the cart position ('my_position') to be our ID in this topology.
        //
        int period = 1;
        MPI_Cart_create( MPI_COMM_WORLD, 1, &worldsz, &period, 1, &cart_comm );
        MPI_Comm_rank( cart_comm, &cart_id );
        MPI_Cart_coords( cart_comm, cart_id, 1, &my_position );
#endif
    }

    /** Destructor frees the allocated one-sided windows.
     */
    ~GameOfLife_Lab3 ()
    {
        MPI_Barrier( MPI_COMM_WORLD );
        MPI_Win_free( &top_win );
        MPI_Win_free( &bottom_win );
    }

protected: //////////////////// VIRTUAL MEMBERS //////////////////////////

    /** Initializes elements of oldc to 0 or 1.
     *  Populate only the local i's on the current rank.
     */
    virtual void InitializeCells ()
    {
        const int offset = my_position * NI;
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

        ////////////////////////////////////////////////////////////////////
        #ifdef TOPOLOGY
            int above_rank, below_rank;
            MPI_Cart_shift( cart_comm, 0, -1, &cart_id, &above_rank );
            MPI_Cart_shift( cart_comm, 0, 1,  &cart_id, &below_rank );
        #else
            int above_rank = myrank <= 0 ? worldsz - 1 : myrank - 1;
            int below_rank = myrank >= worldsz - 1 ? 0 : myrank + 1;
        #endif

        // Use one sided communication to move row from above and
        // below into ghost cells

        // Read row from bottom of above process into top row of ghost cells
        // remember a fence call is needed before and after the get
        //
        MPI_Win_fence( 0, bottom_win );
        MPI_Get( &oldc[0][0], NJ+2, MPI_INT, above_rank, 0, NJ+2, MPI_INT, bottom_win );
        MPI_Win_fence( 0, bottom_win );

        // Read row from below into bottom row of ghost cells
        // remember a fence call is needed before and after the get
        //
        MPI_Win_fence( 0, top_win );
        MPI_Get( &oldc[NI+1][0], NJ+2, MPI_INT, below_rank, 0, NJ+2, MPI_INT, top_win );
        MPI_Win_fence( 0, top_win );
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

    GameOfLife_Lab3( worldsz, myrank, total_NI, total_NJ, debug )
        .Iterate( NSTEPS );

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize ();

    return 0;
}
