/**
    Conway's Game of Life

    Serial implementation and the common MPI stub.
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <mpi.h>

///////////////////////////////////////////////////////////////////////////////
/** Implements the serial version of Conway's Game of Life.
 *  This class can be instantiated on its own or inherited.
 */
class GameOfLife
{
protected:

    int   NI;     //!< Local NI
    int   NJ;     //!< Local NJ
    int** oldc;   //!< Old states
    int** newc;   //!< New states
    int   stepNo; //!< Current step number (while iterating)

public:

    /** Constructs the chunk (allocates memory and initialized arrays).
     */
    GameOfLife
    (
        int ni, int nj
    )
        : NI( ni ), NJ( nj )
        , oldc( NULL ), newc( NULL ), stepNo( 0 )
    {
        if( NI * NJ <= 0 ) {
            return;
        }

        // Allocate the arrays
        //
        oldc = new int*[ NI + 2 ];
        newc = new int*[ NI + 2 ];
        oldc[0] = new int[ ( NI + 2 ) * ( NJ + 2 ) ];
        newc[0] = new int[ ( NI + 2 ) * ( NJ + 2 ) ];

        for( int i = 1; i < NI + 2; ++i )
        {
            oldc[i] = oldc[0] + i * ( NJ + 2 );
            newc[i] = newc[0] + i * ( NJ + 2 );
        }

        // Clear memory
        //
        for( int i = 0; i < NI + 2; ++ i )
        {
            for( int j = 0; j < NJ + 2; ++j )
            {
                oldc[i][j] = 0;
            }
        }
    }

    ~GameOfLife ()
    {
        if( NI * NJ > 0 )
        {
            delete[] oldc[0];
            delete[] newc[0];
            delete[] oldc;
            delete[] newc;
        }
    }

    /**  Play the game of life
     */
    void Iterate( int NSTEPS )
    {
        InitializeCells ();

        for( stepNo = 0; stepNo < NSTEPS; ++stepNo )
        {
            // Fix the  boundary conditions
            //
            exchangeBoundary_LeftRight ();
            exchangeBoundary_TopBottom ();
            /*
                Boundary conditions example:

                NI = 4, NJ = 2           fix LR     mpi TB
                                           ->         ->
                    0 1 2 3        * * * *    * * * *    h G H g
                 0  h g h g        * A B *    b A B a    b A B a
                 1  b A B a        * C D *    d C D c    d C D c
                 2  d C D c        * * * *    * * * *    f E F e
                 3  f E F e
                 4  h G H g        * * * *    * * * *    d C D c
                 5  b a b a        * E F *    f E F e    f E F e
                                   * G H *    h G H g    h G H g
                                   * * * *    * * * *    b A B a
            */

            // Update the cell states
            //
            for( int i = 1; i <= NI; ++i )
            {
                for( int j = 1; j <= NJ; ++j )
                {
                    int im = i - 1;
                    int ip = i + 1;
                    int jm = j - 1;
                    int jp = j + 1;

                    int nsum = oldc[im][jp] + oldc[i][jp] + oldc[ip][jp]
                             + oldc[im][j ]               + oldc[ip][j ]
                             + oldc[im][jm] + oldc[i][jm] + oldc[ip][jm];

                    switch( nsum )
                    {
                        case 3:  newc[i][j] = 1;          break;
                        case 2:  newc[i][j] = oldc[i][j]; break;
                        default: newc[i][j] = 0;
                    }
                }
            }

            // Copy the new state into the old state
            //
            for( int i = 1; i <= NI; ++i )
            {
                for( int j = 0; j <= NJ; ++j )
                {
                    oldc[i][j] = newc[i][j];
                }
            }
        }

        FinalizeEvolution ();
    }

    /** Dump the cell states to an output stream.
     */
    void dumpCells
    (
        std::ostream& of,
        const std::string& prefix = "",
        bool withGhosts = false
    )
    {
        int i1 = 1, i2 = NI, j1 = 1, j2 = NJ;
        if ( withGhosts ) {
            --i1; ++i2; --j1; ++j2;
        }

        for( int i = i1; i <= i2; ++i )
        {
            if ( i == NI + 1 ) {
                of << prefix << "--";
                for( int j = 0; j < NJ + 2; ++j ) {
                    of << "-";
                }
                of << std::endl;
            }

            of << prefix;
            for( int j = j1; j <= j2; ++j )
            {
                of << ( j == NJ + 1 ? "|" : "" );
                of << ( oldc[i][j] ? "O" : " " );
                of << ( j == 0 ? "|" : "" );
            }
            of << std::endl;

            if ( i == 0 ) {
                of << prefix << "--";
                for( int j = 0; j < NJ + 2; ++j ) {
                    of << "-";
                }
                of << std::endl;
            }
        }
    }

    /** Dump the cell states to a file.
     */
    void dumpCells
    (
        const std::string& filename,
        const std::string& prefix = "",
        bool withGhosts = false
    )
    {
        std::cout << "Saved " << filename << std::endl;
        std::ofstream of( filename );
        dumpCells( of );
    }

protected: //////////////////// VIRTUAL MEMBERS //////////////////////////

    /** Initializes elements of oldc to 0 or 1.
     */
    virtual void InitializeCells ()
    {
        // Populate only the local i's on the current rank.
        //
        for( int i = 1; i <= NI; ++i )
        {
            for( int j = 1; j <= NJ; ++j )
            {
                // rand() must be here in the loop to give the same values
                // across all the ranks.
                //
                float x = rand() / ( float(RAND_MAX) + 1);
                if( x >= 0.5 )
                {
                    oldc[i][j] = 1;
                }
            }
        }
    }

    /** Exchanges the left-right boundary of the chunk.
     */
    virtual void exchangeBoundary_LeftRight ()
    {
        // left-right boundary conditions
        //
        for( int i = 0; i <= NI + 1; ++i )
        {
            oldc[i][0]    = oldc[i][NJ];
            oldc[i][NJ+1] = oldc[i][1];
        }
    }

    /** Exchanges the top-bottom boundary of the chunk.
     */
    virtual void exchangeBoundary_TopBottom ()
    {
        // top-bottom boundary conditions
        //
        for( int j = 0; j <= NJ + 1; ++j )
        {
            oldc[0][j]    = oldc[NI][j];
            oldc[NI+1][j] = oldc[1][j];
        }
    }

    /** Iterations are done. Sum the number of live cells and
     *  save the file(s) depicting the obtained life.
     */
    virtual void FinalizeEvolution ()
    {
    }
};

///////////////////////////////////////////////////////////////////////////////
/** Implements the partitioned MPI version of Conway's Game of Life.
 *  This class is not standalone, it should be inherited.
 */
class GameOfLife_MPI : public GameOfLife
{
protected:

    int worldsz; //!< MPI world size
    int myrank;  //!< My rank in the mpi world
    bool debug;  //!< If enabled, writes additional information

public:

    /** Sets up the chunk and stores the mpi world size and our rank.
     */
    GameOfLife_MPI
    (
        int mpiProcs, int mpiRank,
        int ni, int nj, bool debugFlag
    )
        : GameOfLife( ni, nj )
        , worldsz( mpiProcs ), myrank( mpiRank ), debug( debugFlag )
    {
    }

    /** The iterations are done. Sum the number of live cells and
     *  save the file(s) depicting the obtained life.
     */
    void FinalizeEvolution ()
    {
        dumpCells(
            worldsz == 1 ? std::string("life.gol")
                : std::string("life-") + std::to_string(myrank) + ".gol"
        );
        MPI_Barrier( MPI_COMM_WORLD );

        int sum_alive = 0;
        for( int i = 1; i <= NI; ++i )
        {
            for( int j = 1; j <= NJ; ++j )
            {
                sum_alive += newc[i][j];
            }
        }

        // Print final number of live cells.
        //
        if ( worldsz > 1 )
        {
            int local_sum = sum_alive;
            MPI_Reduce( &local_sum, &sum_alive, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
        }

        if ( myrank == 0 ) {
            std::cout << "Number of live cells = " << sum_alive
                << std::endl << std::endl;
        }
    }
};
