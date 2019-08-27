/**
    Calculate pi

    # Lab 2 Exercise 1: Calculate pi using collectives

    Calculates pi using a "dartboard" algorithm. If you're unfamiliar with
    this algorithm, checkout the Wikipedia page on [Monte Carlo
    Integration](http://en.wikipedia.org/wiki/Monte_Carlo_Integration) or
    *Fox et al.(1988) Solving Problems on Concurrent Processors, vol. 1,
    page 207.*

    Hint: All processes should contribute to the calculation, with the
    master averaging the values for pi. Consider using `mpi_reduce` to
    collect results.
*/

#include <iostream>
#include <iomanip>
#include <random>
#include <mpi.h>

static inline double sqr( double x ) { return x * x; }

static std::default_random_engine generator;
static std::uniform_real_distribution<double> distribution(0.0,1.0);

static double dboard( int darts )
{
    // Throw darts at the board
    //
    unsigned long score = 0;
    for( int n = 1; n <= darts; ++n )
    {
        // Generate random numbers for x and y coordinates
        //
        double x_coord = 2.0 * distribution(generator) - 1.0;
        double y_coord = 2.0 * distribution(generator) - 1.0;

        // If dart lands in circle, increment score
        //
        if( sqr(x_coord) + sqr(y_coord) <= 1.0 ) {
            ++score;
        }
    }

    return 4.0 * double(score) / darts;
}

int main( int argc, char** argv )
{
    /** Number of throws at dartboard  */
    int DARTS = argc >= 2 ? atoi( argv[1] ) : 50000;

    /** Number of times "darts" is iterated */
    int ROUNDS = argc >= 10 ? atoi( argv[1] ) : 10;

    // Setup MPI
    //
    MPI_Init( &argc, &argv );

    int myrank, worldsz;
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    MPI_Comm_size( MPI_COMM_WORLD, &worldsz );

    std::cout << "MPI task " << myrank << " has started..." << std::endl;
    if( myrank == 0 ) {
        std::cout << "Using " << worldsz << " tasks to compute pi" << std::endl;
    }

    generator.seed( myrank );

    double avg_pi = 0, pi_sum = 0;
    for( int i = 0; i < ROUNDS; ++i )
    {
        double my_pi = dboard( DARTS );

        // Use MPI_SUM reduction function.
        //
        int rc = MPI_Reduce( &my_pi, &pi_sum, 1, MPI_DOUBLE, MPI_SUM,
                             0, MPI_COMM_WORLD);
        if( rc != MPI_SUCCESS ) {
            std::cerr << myrank << ": Reduce failed, round = " << i << std::endl;
        }

        if( myrank == 0 ) // Master will printout the result
        {
            // Update the average pi
            //
            double pi = pi_sum / worldsz;
            avg_pi = ( avg_pi * i + pi ) / ( i + 1 );

            std::cout << DARTS * ( i + 1 ) << " throws, average pi = "
                << std::setprecision(14) << avg_pi << std::endl;
        }
    }

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize ();

    return 0;
}
