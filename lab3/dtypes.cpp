/**
    # Exercise 1 - MPI I/O

    MPI I/O is used so that results can be written to the same file in parallel.
    Take the serial hello world programs and modify them so that instead of
    writing the output to screen the output is written to a file using MPI I/O.

    The simplest solution is likely to be for you to create a character buffer,
    and then use the MPI_File_write_at function.

    STL file format:
    ```
        UINT8[80] – Header
        UINT32 – Number of triangles

        foreach triangle
        REAL32[3] – Normal vector
        REAL32[3] – Vertex 1
        REAL32[3] – Vertex 2
        REAL32[3] – Vertex 3
        UINT16 – Attribute byte count
        end
    ```
    @see: https://en.wikipedia.org/wiki/STL_(file_format
*/

#include <iostream>
#include <iomanip>

#include <cstdio>
#include <cstdint>
#include <cstring>

#include <mpi.h>



class STL_Model
{
public:

    enum { HDR_SIZE = 80 };

protected:

    struct __attribute__((packed)) Triangle
    {
        float n[3];         //!< Normal vector
        float v1[3];        //!< Vertex 1
        float v2[3];        //!< Vertex 2
        float v3[3];        //!< Vertex 3 */
        uint16_t attrib;    //!< Attribute byte count
    };

    uint32_t   n_tri;           //!< Number of triangles
    Triangle*  tri;             //!< Triangles
    char       hdr[HDR_SIZE];   //!< Header

public:

    STL_Model ()
        : n_tri( 0 ), tri( NULL )
    {
    }

    ~STL_Model ()
    {
        if( tri != NULL ) {
            delete[] tri;
        }
    }

    virtual void read( const char* fname )
    {
        int pe_size, pe_rank;
        MPI_Comm_size( MPI_COMM_WORLD, &pe_size );
        MPI_Comm_rank( MPI_COMM_WORLD, &pe_rank );

        FILE* fp = fopen( fname, "rb" );

        if( pe_rank == 0 ) {
            std::cout << std::endl 
                << "Reading STL file: " << fname
                << std::endl << "Header size =  " << HDR_SIZE 
                << std::endl << "sizeof(Triangle) = " << sizeof(Triangle )
                << std::endl;
        }

        // Read STL header
        //
        size_t rc = fread( hdr, 1, HDR_SIZE, fp);

        // Make sure it's a binary STL file with the correct header length
        //
        if ( rc != HDR_SIZE || strncmp( hdr, "solid", 5) == 0)
        {
            std:: cerr << "ASCII STL files are not supported!" << std::endl;
            exit( -1 );
        }

        // Read how m   any triangles the file contains
        //
        rc = fread( &n_tri, sizeof(uint32_t), 1, fp);
        if( pe_rank == 0 ) {
            std::cout << "Found: " << n_tri << " triangles" << std::endl;
        }

        // Allocate memory for triangles, and read them
        //
        tri = new Triangle[ n_tri ];
        rc = fread( tri, sizeof(Triangle), n_tri, fp );
        if( rc != n_tri ) {
            std:: cerr << "The STL file is too short!" << std::endl;
            exit( -1 );
        }

        fclose( fp );
        if (pe_rank == 0) {
            std::cout << "Done reading." << std::endl;
        }
    }

    virtual void write( const char* fname )
    {
        int pe_size, pe_rank;
        MPI_Comm_size( MPI_COMM_WORLD, &pe_size );
        MPI_Comm_rank( MPI_COMM_WORLD, &pe_rank );

        // Only the master process can write the complete file
        //
        if ( pe_rank != 0 ) 
        {  
            return;
        }

        FILE* fp = fopen( fname, "wb" );
        if( pe_rank == 0 ) {
            std::cout << std::endl << "Writing STL file: " << fname << std::endl;
        }

        // Write STL header
        //
        fwrite( hdr, 1, HDR_SIZE, fp );

        // Write number of triangles
        //
        fwrite( &n_tri, sizeof(uint32_t), 1, fp );

        // Write all triangles
        //
        fwrite( tri, sizeof(Triangle), n_tri, fp );

        fclose( fp );

        std::cout << "Done writing." << std::endl;
    }
};

class STL_Model_MPI : public STL_Model
{
    void debug_write( int pe_rank, MPI_Offset offset, size_t len )
    {
        std::cout << std::setw(2) << pe_rank << ": written "  \
            << std::setw(7) << offset << " -- "  \
            << std::setw(7) << offset + len - 1 \
            << ", len = " << len << " octets" << std::endl;
    }

public:

    virtual void write( const char* fname )
    {
        int pe_size, pe_rank;
        MPI_Comm_size( MPI_COMM_WORLD, &pe_size );
        MPI_Comm_rank( MPI_COMM_WORLD, &pe_rank );

        // Open the file for writing
        //
        MPI_File fp;
        MPI_File_open( MPI_COMM_WORLD,
                   fname, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                   MPI_INFO_NULL, &fp );

        if( pe_rank == 0 ) {
            std::cout << std::endl << "MPI Writing STL file: " << fname << std::endl;
        }

        // Write STL header
        //
        // fwrite( hdr, 1, HDR_SIZE, fp );
        //
        MPI_Offset offset = 0;
        size_t len = HDR_SIZE;
        if( pe_rank == 0 ) 
        {
            MPI_Status status;
            MPI_File_write_at( fp, offset, hdr, len, MPI_CHAR, &status );
            debug_write( pe_rank, offset, len );
        }
        offset += len;

        // Write number of triangles
        //
        // fwrite( &n_tri, sizeof(uint32_t), 1, fp );
        //
        len = sizeof(uint32_t);
        if( pe_rank == 0 ) 
        {
            MPI_Status status;
            MPI_File_write_at( fp, offset, &n_tri, len, MPI_CHAR, &status );
            debug_write( pe_rank, offset, len );
        }
        offset += len;

        MPI_Barrier( MPI_COMM_WORLD ); // Wait for the header to be written

        // Write all triangles
        //
        // fwrite( tri, sizeof(Triangle), n_tri, fp );
        //
        size_t chunk_size = n_tri / pe_size;
        offset += pe_rank * chunk_size * sizeof(Triangle);
        len = ( pe_rank < pe_size - 1 ? chunk_size : n_tri - pe_rank * chunk_size ) 
            * sizeof(Triangle);

        MPI_Status status;
        MPI_File_write_at_all( fp, offset, &tri[ pe_rank * chunk_size ], len, MPI_CHAR, &status );
        debug_write( pe_rank, offset, len );

        // Close the file
        //
        MPI_File_close( &fp );

        if( pe_rank == 0 ) {
            std::cout << "Done writing." << std::endl;
        }
    }
};

int main( int argc, char** argv )
{
    MPI_Init( &argc, &argv );

    STL_Model_MPI model;

    model.read ( argc >= 2 ? argv[1] : "data/sphere.stl" );
    model.write( argc >= 3 ? argv[2] : "out.stl"         );

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize ();

    return 0;
}


