/**
    Calculate pi

    # Lab 2 Exercise 3: Find pi using non-blocking communications

    Use a non-blocking send to try to overlap communication and computation.
    Take the code from Exercise 1 as your starting point.
*/

#include <iostream>
#include <iomanip>
#include <random>
#include <mpi.h>

static inline double sqr( double x ) { return x * x; }

static double dboard( int seed, int darts )
{
    std::default_random_engine generator( seed );
    std::uniform_real_distribution<double> distribution(0.0,1.0);

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

    double avg_pi = 0;
    for( int i = 0; i < ROUNDS; ++i )
    {
        double my_pi = dboard( i, DARTS );

        if( myrank != 0 )  // Slave
        {
            MPI_Request request;
            int rc = MPI_Isend( &my_pi, 1, MPI_DOUBLE,
                                0, i, MPI_COMM_WORLD, &request );

            if( rc != MPI_SUCCESS ) {
                std::cerr << myrank << ": ISend failed; round = " << i << std::endl;
            }

            if( i >= ROUNDS - 1 ) // Wait at the last round
            {
                MPI_Status status;
                MPI_Wait( &request, &status );
            }
        }
        else // Master
        {
            // Get sum ov pi from the slaves
            //
            double slave_pi_sum = 0;
            for( int n = 1; n < worldsz; ++n )
            {
                MPI_Status status;
                double recv_pi = 0;
                int rc = MPI_Recv( &recv_pi, 1, MPI_DOUBLE, MPI_ANY_SOURCE,
                                   i, MPI_COMM_WORLD, &status);
                if( rc != MPI_SUCCESS ) {
                    std::cerr << myrank << ": Receive failed, round = " << i << std::endl;
                }
                slave_pi_sum += recv_pi;
            }

            // Update the average pi
            //
            double pi = ( slave_pi_sum + my_pi ) / worldsz;
            avg_pi = ( avg_pi * i + pi ) / ( i + 1 );

            std::cout << DARTS * ( i + 1 ) << " throws, average pi = "
                << std::setprecision(14) << avg_pi << std::endl;
        }
    }

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize ();

    return 0;
}
