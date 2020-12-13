#include <stdio.h>
#include <mpi.h>
#include <malloc.h>
#include <stdbool.h>
#include <stdlib.h>


void kernel (double * K , int dim_kernel , int type_kernel , int weight) ;

int main (int argc , char * argv[]){
    
    //Initialize matrix for the image and the kernel
    if (argc <=1){
        fprintf (stderr , " Required at least 3 parameters : PGM Image ,"
         "type of the kernel , dimension of the kernel and weight (if needed) ") ;
        
        exit (-1) ;
    }
    int N[2] = {256,256} ;
    char * image = argv[1] ;
    short int * M ;
    M = (short int*)malloc(N[0] * N[1] * sizeof(short int));
    FILE *f = fopen (image , "rb") ;
    fread (M , 1 , N[0] * N[1] , f) ;
    fclose (f) ;



    int type_kernel ;
    int dim_kernel ;
    type_kernel = atoi (argv[2]) ;
    dim_kernel = atoi (argv[3]) ;
    int weight ;
    if (atoi (argv[4])) { weight = atoi (argv[4]) ;}
    else {weight = 0 ;}

    double * K ;
    K = (double*)malloc( dim_kernel * dim_kernel * sizeof(double)) ;
    kernel ( K , dim_kernel , type_kernel , weight) ;

    
    MPI_Init(&argc , &argv);
    MPI_Status status;

    //Create a new comm world
    int size;
    MPI_Comm_size(MPI_COMM_WORLD , &size);

    int dims[2] = {0,0};
    MPI_Dims_create(size ,2 , dims);

    int periods[2] = {false , false};
    int reorder = true;

    MPI_Comm new;
    MPI_Cart_create(MPI_COMM_WORLD ,2 , dims ,
    periods , reorder , &new);

    int rank;
    MPI_Comm_rank(new , &rank);



    if (rank == 0){
        for (int i = 0 ; i < dim_kernel ; i++){
            for (int j = 0 ; j < dim_kernel ; j++)
                printf (" %9.2f ",K[i * dim_kernel + j]) ;
            printf ("\n") ;
        }    
    }

    //Define the number of elements in each dimension for each proc
    int local_dim[2];
    
    for( int i = 0; i < 2 ; i++ ){ 

    local_dim[i] = N[i]/dims[i];
    if( rank == (size-1)  )    local_dim[i] += N[i]%dims[i];

    }

    double * local_M ;
    local_M = (double*)malloc(local_dim[0] * local_dim[1] * sizeof(double)) ;

    //Initialize the local matrix
    /*int offset = rank * local_dim[0] ;
    for (int i = 0 ; i < local_dim[0] ; i++)
        for(int j = 0 ; j < local_dim[1] ; j++)
            local_M[i * local_dim[0] + j] = M[(offset + i) * N[0] + (j + offset)]*/

    free ( local_M );
    MPI_Finalize () ;

    free ( M ) ;
    free( K ) ;

    return 0 ;

}



void kernel (double * K, int dim_kernel , int type_kernel , int weight){
    //Mean Kernel = 0
    //Weighted Kernel = 1
    //Guassian Kernel = 2
    switch (type_kernel)
    {
    case 0 : 

        for(int i = 0 ; i < dim_kernel ; i++)
            for(int j = 0 ; j < dim_kernel ; j++)
                K[i*dim_kernel + j] = 1 / (float)(dim_kernel * dim_kernel);
        
        break;

    case 1 :

        for(int i = 0 ; i < dim_kernel ; i++)
            for(int j = 0 ; j < dim_kernel ; j++){
                if(i == ((dim_kernel-1) / 2) && j == i)
                    K[i*dim_kernel + j] = weight ;
                else
                    K[i*dim_kernel + j] = (1 - weight)/(dim_kernel^2 -1) ;
            }

        break;
                
    case 2 :

        for(int i = 0 ; i < dim_kernel ; i++)
            for(int j = 0 ; j < dim_kernel ; j++)
                K[i*dim_kernel + j] = 0.2 ;
    
    default:
        break;
    }

    
}