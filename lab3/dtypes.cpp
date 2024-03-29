/**
    Exercise 2 - MPI I/O and derived types

    Take the serial stl reader and modify it such that the stl file is read 
    (and written) in parallel using collective MPI I/O. Use derived types such 
    that the file can be read/written with a maximum of 3 I/O operations per 
    read and write.

    The simplest solution is likely to create a derived type for each triangle,
    and then use the MPI_File_XXXX_at_all function. A correct solution will 
    have the same MD5 hash for both stl models (input and output), unless 
    the order of the triangles has been changed.

    md5sum out.stl data/sphere.stl

        822aba6dc20cc0421f92ad50df95464c  out.stl
        822aba6dc20cc0421f92ad50df95464c  data/sphere.stl


    STL file format:    (all fields are little endian in the file)
    ```
        UINT8[80]       – Header
        UINT32          – Number of triangles

        foreach triangle
            REAL32[3]   – Normal vector
            REAL32[3]   – Vertex 1
            REAL32[3]   – Vertex 2
            REAL32[3]   – Vertex 3
            UINT16      – Attribute byte count
        end
    ```
    @see: https://en.wikipedia.org/wiki/STL_(file_format

    @warning: This code works only on a little endian machine.
*/

#include <iostream>
#include <iomanip>

#include <cstdio>
#include <cstdint>
#include <cstring>

#include <mpi.h>

///////////////////////////////////////////////////////////////////////////////

class STL_Model
{
public:

    enum { HDR_SIZE = 80 };

protected:

    /** A triangle STL element
     */
    struct __attribute__((packed)) Triangle
    {
        float    n  [3];    //!< Normal vector
        float    v1 [3];    //!< Vertex 1
        float    v2 [3];    //!< Vertex 2
        float    v3 [3];    //!< Vertex 3 */
        uint16_t attrib;    //!< Attribute byte count
    };

    char      hdr[HDR_SIZE];   //!< Header
    uint32_t  n_tri;           //!< Number of triangles

    Triangle* tri;             //!< Triangles

    ///////////////////////////////////////////////////////////////////////////

    bool debugMode;

    void debug( const char* info, MPI::Offset offset, size_t len )
    {
        if( debugMode )
        {
            std::cout << std::setw(2) << MPI::COMM_WORLD.Get_rank () 
                << ": " << info << " "  \
                << std::setw(7) << offset << " -- "  \
                << std::setw(7) << offset + len - 1 \
                << ", len = " << len << " octets" << std::endl;
        }
    }

public:

    STL_Model ()
        : n_tri( 0 ), tri( NULL ), debugMode( false )
    {
    }

    ~STL_Model ()
    {
        if( tri != NULL ) {
            delete[] tri;
        }
    }

    void setDebug( bool flag )
    {
        debugMode = flag;
    }

    virtual void read( const char* fname )
    {
        FILE* fp = fopen( fname, "rb" );

        int pe_rank = MPI::COMM_WORLD.Get_rank ();
        if( pe_rank == 0 && debugMode ) {
            std::cout << std::endl 
                << "Reading STL file: " << fname << std::endl;
            std::cout << std::endl << "STL header size = " << HDR_SIZE 
                << std::endl << "sizeof(Triangle) = " << sizeof(Triangle )
                << std::endl;
        }

        // Read STL header
        //
        size_t rc = fread( hdr, 1, HDR_SIZE, fp);

        // Make sure it's a binary STL file with the correct header length
        //
        if ( rc != HDR_SIZE || strncmp( hdr, "solid", 5 ) == 0)
        {
            std:: cerr << "ASCII STL files are not supported!" << std::endl;
            exit( -1 );
        }

        // Read how m   any triangles the file contains
        //
        rc = fread( &n_tri, sizeof(uint32_t), 1, fp);
        if( pe_rank == 0 && debugMode ) {
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
        if( pe_rank == 0 && debugMode ) {
            std::cout << "Done reading." << std::endl;
        }
    }

    virtual void write( const char* fname )
    {
        // Only the master process can write the complete file
        //
        int pe_rank = MPI::COMM_WORLD.Get_rank ();
        if ( pe_rank != 0 ) 
        {  
            return;
        }

        FILE* fp = fopen( fname, "wb" );
        if( pe_rank == 0 && debugMode ) {
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

        if( pe_rank == 0 && debugMode ) {
            std::cout << "Done writing." << std::endl;
        }
    }
};

///////////////////////////////////////////////////////////////////////////////

/** An implementation where write is using MPI::BYTE
 */
class STL_Model_MpiByte : public STL_Model
{
public:

    virtual void write( const char* fname )
    {
        // Open the file for writing
        //
        MPI::File fp = MPI::File::Open( MPI::COMM_WORLD,
                   fname, MPI::MODE_CREATE | MPI::MODE_WRONLY,
                   MPI::INFO_NULL );

        int pe_rank = MPI::COMM_WORLD.Get_rank ();
        if( pe_rank == 0 && debugMode ) {
            std::cout << std::endl << "MPI (Simple) Writing STL file: " 
                << fname << std::endl;
        }

        // Write STL header
        //
        MPI::Offset offset = 0;
        size_t len = HDR_SIZE;
        if( pe_rank == 0 ) 
        {
            MPI::Status status;
            fp.Write_at( offset, hdr, len, MPI::BYTE, status );
            debug( "write", offset, len );
        }
        offset += len;

        // Write number of triangles
        //
        len = sizeof(uint32_t);
        if( pe_rank == 0 ) 
        {
            MPI::Status status;
            fp.Write_at( offset, &n_tri, len, MPI::BYTE, status );
            debug( "write", offset, len );
        }
        offset += len;

        MPI::COMM_WORLD.Barrier (); // Wait for the header to be written

        // Write all triangles
        //
        int pe_size = MPI::COMM_WORLD.Get_size ();
        size_t chunk_size = n_tri / pe_size;
        offset += pe_rank * chunk_size * sizeof(Triangle);
        len = ( pe_rank < pe_size - 1 ? chunk_size : n_tri - pe_rank * chunk_size ) 
            * sizeof(Triangle);

        MPI::Status status;
        fp.Write_at_all( offset, &tri[ pe_rank * chunk_size ], len, MPI::BYTE, status );
        debug( "write", offset, len );

        // Close the file
        //
        fp.Close ();

        if( pe_rank == 0 && debugMode ) {
            std::cout << "Done writing." << std::endl;
        }
    }
};

///////////////////////////////////////////////////////////////////////////////

/** An implementation where both read and write are using MPI:Datatype
 */
class STL_Model_MpiDT : public STL_Model
{
    MPI::Datatype MPI_Triangle;
    MPI::Datatype MPI_TrianglePacked;

public:

    STL_Model_MpiDT ()
    {
        int blockcounts[2];
        MPI::Aint offsets[2];
        MPI::Datatype oldtypes[2];

        // Need to first figure offset by getting size of MPI_FLOAT
        //
        MPI::Aint lb, extent;
        MPI::FLOAT.Get_extent( lb, extent );

        // Setup description of the 4*3 MPI_FLOAT
        // This block covers the fields: n[3], v1[3], v2[3], and v3[3]
        //
        offsets     [0] = 0;
        oldtypes    [0] = MPI::FLOAT;
        blockcounts [0] = 4 * 3;

        // Setup description of the 1 MPI_INT field
        // This block covers the field `attrib`
        //
        offsets     [1] = 4 * 3 * extent;
        oldtypes    [1] = MPI::UNSIGNED_SHORT;
        blockcounts [1] = 1;

        // Define structured type and commit it
        //
        MPI_Triangle = MPI::Datatype::Create_struct( 2, blockcounts, offsets, oldtypes );
        MPI_Triangle.Commit ();

        // Since wi have packed structure, we need to reduce the extent
        //
        MPI_Triangle.Get_extent( lb, extent );
        MPI_TrianglePacked = MPI_Triangle.Create_resized( lb, extent -2 );
        MPI_TrianglePacked.Commit ();

        // Report creation
        //
        if ( MPI::COMM_WORLD.Get_rank() == 0 && debugMode ) {
            std::cout << "MPI_Triangle extent = " << extent << std::endl;
            MPI_TrianglePacked.Get_extent( lb, extent );
            std::cout << "MPI_TrianglePacked extent = " << extent << std::endl;
        }
    }

    virtual void read( const char* fname )
    {
        // Open the file for reading
        //
        MPI::File fp = MPI::File::Open( MPI::COMM_WORLD,
                   fname, MPI::MODE_RDONLY,
                   MPI::INFO_NULL );

        int pe_rank = MPI::COMM_WORLD.Get_rank ();
        if( pe_rank == 0 && debugMode ) {
            std::cout << std::endl 
                << "MPI (DT) Reading STL file: " << fname << std::endl;
            std::cout << std::endl << "STL header size = " << HDR_SIZE 
                << std::endl << "sizeof(Triangle) = " << sizeof(Triangle )
                << std::endl;
        }

        // Read STL header
        //
        MPI::Offset offset = 0;
        size_t len = HDR_SIZE;
        // if( pe_rank == 0 ) // everyone must read header
        {
            MPI::Status status;
            fp.Read_at_all( offset, hdr, len, MPI::BYTE, status );
            debug( "readH", offset, len );

            // Make sure it's a binary STL file with the correct header length
            //
            if ( strncmp( hdr, "solid", 5 ) == 0)
            {
                std:: cerr << "ASCII STL files are not supported!" << std::endl;
                exit( -1 );
            }
        }
        offset += len;

        // Read number of triangles
        //
        len = sizeof(uint32_t);
        // if( pe_rank == 0 ) // everyone must read header
        {
            MPI::Status status;
            fp.Read_at_all( offset, &n_tri, 1, MPI::INT, status );
            debug( "readN", offset, len );

            // Allocate memory for triangles
            //
            tri = new Triangle[ n_tri ];
        }
        offset += len;

        if( pe_rank == 0 && debugMode ) {
            std::cout << "Found: " << n_tri << " triangles" << std::endl;
        }

        // Read all triangles using MPI_TrianglePacked datatype
        //
        int pe_size = MPI::COMM_WORLD.Get_size ();
        size_t chunk_size = n_tri / pe_size;
        offset += pe_rank * chunk_size * sizeof(Triangle);
        len = ( pe_rank < pe_size - 1 ? chunk_size : n_tri - pe_rank * chunk_size );

        MPI::Status status;
        fp.Read_at_all( offset, &tri[ pe_rank * chunk_size ], len, 
                         MPI_TrianglePacked, status );
        debug( "readD", offset, len * sizeof(Triangle) );

        // Close the file
        //
        fp.Close ();

        if( pe_rank == 0 && debugMode ) {
            std::cout << "Done reading." << std::endl;
        }

        // MPI::COMM_WORLD.Barrier (); // Wait for the header to be read
    }

    virtual void write( const char* fname )
    {
        // Open the file for writing
        //
        MPI::File fp = MPI::File::Open( MPI::COMM_WORLD,
                   fname, MPI::MODE_CREATE | MPI::MODE_WRONLY,
                   MPI::INFO_NULL );

        int pe_rank = MPI::COMM_WORLD.Get_rank ();
        if( pe_rank == 0 && debugMode ) {
            std::cout << std::endl << "MPI (DT) Writing STL file: "
                << fname << std::endl;
        }

        // Write STL header
        //
        MPI::Offset offset = 0;
        size_t len = HDR_SIZE;
        if( pe_rank == 0 ) 
        {
            MPI::Status status;
            fp.Write_at( offset, hdr, len, MPI::BYTE, status );
            debug( "write", offset, len );
        }
        offset += len;

        // Write number of triangles
        //
        len = sizeof(uint32_t);
        if( pe_rank == 0 ) 
        {
            MPI::Status status;
            fp.Write_at( offset, &n_tri, 1, MPI::INT, status );
            debug( "write", offset, len );
        }
        offset += len;

        MPI::COMM_WORLD.Barrier (); // Wait for the header to be written

        // Write all triangles using MPI_TrianglePacked datatype
        //
        int pe_size = MPI::COMM_WORLD.Get_size ();
        size_t chunk_size = n_tri / pe_size;
        offset += pe_rank * chunk_size * sizeof(Triangle);
        len = ( pe_rank < pe_size - 1 ? chunk_size : n_tri - pe_rank * chunk_size );

        MPI::Status status;
        fp.Write_at_all( offset, &tri[ pe_rank * chunk_size ], len, 
                         MPI_TrianglePacked, status );
        debug( "write", offset, len * sizeof(Triangle) );

        // Close the file
        //
        fp.Close ();

        if( pe_rank == 0 && debugMode ) {
            std::cout << "Done writing." << std::endl;
        }
    }
};


int main( int argc, char** argv )
{
    const char* infname  = argc >= 2 ? argv[1]         : "data/sphere.stl" ;
    const char* outfname = argc >= 3 ? argv[2]         : "out.stl"         ;
    const int   useClass = argc >= 4 ? atoi( argv[3] ) : 2                 ;
    const int   reps     = argc >= 5 ? atoi( argv[4] ) : 1                 ;

    MPI::Init( argc, argv );

    if( MPI::COMM_WORLD.Get_rank () == 0 )
    { 
        std::cout << std::endl << "*** dtypes is running; parameters:" << std::endl
            << "    infname = " << infname << "," << std::endl
            << "    outfname = " << outfname << "," << std::endl
            << "    useClass = " 
            << ( useClass == 2 ? "MPI (DT)" : useClass == 1 ? "MPI (Byte)" : "Serial" )
            << ", reps = " << reps << "," << std::endl
            << "    MPI world size = " << MPI::COMM_WORLD.Get_size () << std::endl;
   }

    double sumT = 0;
    for( int n = 0; n < reps; ++n )
    {
        double T1, T2;
        if( useClass == 2 ) // Test STL_Model_MpiDT
        {
            STL_Model_MpiDT model;
            model.read( infname );
            model.setDebug( reps == 1 );
            T1 = MPI_Wtime ();
            model.write( outfname );
            T2 = MPI_Wtime ();
        }
        else if ( useClass == 1 ) // Test STL_Model_MpiByte
        {
            STL_Model_MpiByte model;
            model.read( infname );
            model.setDebug( reps == 1 );
            T1 = MPI_Wtime ();
            model.write( outfname );
            T2 = MPI_Wtime ();
        }
        else // Test STL_Model (base class, serial)
        {
            STL_Model model;
            model.read( infname );
            model.setDebug( reps == 1 );
            T1 = MPI_Wtime ();
            model.write( outfname );
            T2 = MPI_Wtime ();
        }

        // Calculate round trip time and print
        //
        if( MPI::COMM_WORLD.Get_rank () == 0 )
        {
            double deltaT = T2 - T1;
            sumT += deltaT;

            if( n == 0 ) {
                std::cout << std::endl
                    << "Rep#        T1              T2            deltaT" << std::endl;
            }

            std::cout << std::setw(4) << n << "  " << std::fixed //<< std::setfill( '0' )
                << std::setw(14) << std::setprecision(8) << T1 << "  "
                << std::setw(14) << std::setprecision(8) << T2 << "  "
                << std::setw(2) << std::setprecision(8) << deltaT
                << std::endl;
        }
    }

    if( MPI::COMM_WORLD.Get_rank () == 0 ) 
    {
        double avgT = sumT * 1e6 / reps;
        std::cout << std::endl << std::setprecision(2) << "*** avgT = " << avgT << " us" << std::endl;
    }

    MPI::COMM_WORLD.Barrier ();
    MPI::Finalize ();

    return 0;
}


