#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

int main( int argc , char* argv[] ){

    MPI_Init( &argc , &argv );
    MPI_Status status;


    int size;
    MPI_Comm_size( MPI_COMM_WORLD , &size );

    int dims[2] = {0,0};
    MPI_Dims_create( size ,2 , dims );

    int periods[2] = { false , false };
    int reorder = true;

    MPI_Comm new;
    MPI_Cart_create( MPI_COMM_WORLD ,2 , dims ,
    periods , reorder , &new );

    int rank;
    MPI_Comm_rank( new , &rank );

    
    int N[2];
    //Ask for the dimension of the problem
    if( rank == 0 ){
        
        printf( "Insert problem size on first axis:\n" );
        scanf( "%d" , &N[0] );
        printf( "Insert problem size on second axis:\n" );
        scanf( "%d" , &N[1] );

        
    }

    MPI_Barrier( new );
    MPI_Bcast( N ,2 ,MPI_INT ,0 , new );
    
    
    
    
    //Let's find the neighbours for each proc
    

    //Collect the coordinates
    int my_coord[2];

    MPI_Cart_coords( new , rank , 2 , my_coord );
    

    //Define the number of elements in each dimension for each proc
    int local_dim[2];
    
    for( int i = 0; i < 2 ; i++ ){ 

    local_dim[i] = N[i]/dims[i];
    if( rank == (size-1)  )    local_dim[i] += N[i]%dims[i];

    }
    
    int iEnd = local_dim[0] + 1;
    int jEnd = local_dim[1] + 1;
    
    printf( "[Rank: %d] my dimensions are (%d,%d)\n" , rank , local_dim[0] , local_dim[1] );

    double phi[iEnd][jEnd][2];

    int MaxBufLen = max( local_dim[0] , local_dim[1] );

    double fieldSend[MaxBufLen];
    double fieldRecv[MaxBufLen];

    int t0 = 0;
    int t1 = 1;
    double eps = 10^(-14);
    double maxdelta = 10;
    int iter;
    int maxiter = 10000;
    do while( maxdelta > eps && iter < maxiter ){

        for( int dir = 0 ; dir < 2 ; dir++ ){
            int source;
            int dest;

            MPI_Cart_shift( new , dir , 1 , &source , &dest );

            if( source ){
                MPI_Recv( fieldRecv , local_dim[dir] , MPI_DOUBLE , source , 0 , new , &status );
            }

            if ( dest ){
                
            }
        }


    }

    MPI_Finalize();
    return 0;
}