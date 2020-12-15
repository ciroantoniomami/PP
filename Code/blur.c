#include <stdio.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdlib.h>
#define MPI_Proc_null -1 

void kernel (double * K , int dim_kernel , int type_kernel , int weight) ;

int main (int argc , char * argv[]){
    
    //Initialize matrix for the image and the kernel
    if (argc <=1){
        fprintf (stderr , " Required at least 3 parameters : PGM Image ,"
         "type of the kernel , dimension of the kernel and weight (if needed) ") ;
        
        exit (-1) ;
    }
    int N[2] = {256,256} ;
    short int * M ;
    char * image = argv[1] ;
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
    int kdim = dim_kernel / 2 ;
    double * K ;
    K = (double*)malloc( dim_kernel * dim_kernel * sizeof(double)) ;
    kernel ( K , dim_kernel , type_kernel , weight) ;

    
    MPI_Init(&argc , &argv);
    MPI_Status status[4] ;
    MPI_Request request[4] ;
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

    //Define the number of elements in each dimension for each proc
    int local_dim[2];
    
    local_dim[0] = N[0] / dims[0] ;
    local_dim[1] = N[1] / dims[1] ;

    short int * local_M ;
    local_M = (short int *)calloc(local_dim[0] * local_dim[1] , sizeof(short int)) ;

    MPI_Datatype type, resizedtype ;
    int starts [2] = {0 , 0} ;

    MPI_Type_create_subarray (2, N , local_dim , starts , MPI_ORDER_C , MPI_SHORT , &type) ;
    MPI_Type_create_resized (type, 0, local_dim[0] * sizeof(short int), &resizedtype) ;
    MPI_Type_commit (&resizedtype) ;

    int counts[size] ;
    int displs[size] ;





    if (rank == 0){


        for (int i = 0 ; i < dim_kernel ; i++){
            for (int j = 0 ; j < dim_kernel ; j++)
                printf (" %9.2f ",K[i * dim_kernel + j]) ;
            printf ("\n") ;
        }

    for (int i = 0; i < size; i++) counts[i] = 1 ;
    for (int i = 0; i < size; i++) displs[i] = i % dims[0] + i / dims[0] * local_dim[1] * dims[0] ;

    }
    MPI_Bcast (counts , size , MPI_INT , 0 , new) ;
    MPI_Bcast (displs , size , MPI_INT , 0 , new) ;

    MPI_Barrier (new) ;
    /* send submatrices to all processes */
    MPI_Scatterv (M , counts , displs , resizedtype , local_M ,
     local_dim[0] * local_dim[1] , MPI_SHORT, 0, new) ;

    MPI_Barrier (new) ;
   /* for (int k = 0 ; k < size ; k++){
        if ( rank == k){
    for (int i = 0 ; i < local_dim[1] ; i++){
     for (int j = 0 ; j < local_dim[0] ; j++)
         printf (" %d ",local_M[i * local_dim[0] + j]) ;
     printf ("\n") ;
    }
        }
    }*/

    // The structure that will be used to exchange data between   
    short int * ubuffer = (short int *) calloc ( kdim * local_dim[0] , sizeof ( short int)) ;    
    short int * dbuffer = (short int *) calloc ( kdim * local_dim[0] , sizeof ( short int)) ;
    short int * lbuffer = (short int *) calloc ( local_dim[1] * kdim , sizeof ( short int)) ; 
    short int * rbuffer = (short int *) calloc ( local_dim[1] * kdim , sizeof ( short int)) ; 
    
    int subsizes_col[2] = {kdim , local_dim[1]} ;
    int subsizes_row[2] = {local_dim[0] , kdim} ;

    MPI_Datatype type_col , colhalo, type_row , rowhalo ;

    MPI_Type_create_subarray (2, local_dim , subsizes_col , starts , MPI_ORDER_C , MPI_SHORT , &type_col) ;
    MPI_Type_create_resized (type_col, 0, kdim * sizeof(short int), &colhalo) ;
    MPI_Type_commit (&colhalo) ;

    MPI_Type_create_subarray (2, local_dim , subsizes_row , starts , MPI_ORDER_C , MPI_SHORT , &type_row) ;
    MPI_Type_create_resized (type_row, 0, local_dim[0] * sizeof(short int), &rowhalo) ;
    MPI_Type_commit (&rowhalo) ;

    int sourceleft , destright ;

    MPI_Cart_shift (new , 0 , 1 , &sourceleft , &destright) ;

    if (sourceleft != MPI_PROC_NULL){
        MPI_Irecv (lbuffer , 1 , colhalo , sourceleft , 123 , new , &request[0]) ;
    }

    if (destright != MPI_PROC_NULL){
        MPI_Send (&local_M[ local_dim[0] - kdim  ] , 1 , colhalo , destright , 123 , new) ;
    }

    if(sourceleft != MPI_PROC_NULL){
        MPI_Wait ( &request[0] , &status[0]) ;
    }

    int sourceright , destleft ;

    MPI_Cart_shift (new , 0 , -1 , &sourceright , &destleft) ;

    if (sourceright != MPI_PROC_NULL){
        MPI_Irecv (rbuffer , 1 , colhalo , sourceright , 124 , new , &request[1]) ;
    }

    if (destleft != MPI_PROC_NULL){
        MPI_Send (local_M , 1 , colhalo , destleft , 124 , new) ;
    }

    if(sourceright != MPI_PROC_NULL){
        MPI_Wait ( &request[1] , &status[1]) ;
    }

    int sourceup , destdown ;

    MPI_Cart_shift (new , 1 , 1 , &sourceup , &destdown) ;

    if (sourceup != MPI_PROC_NULL){
        MPI_Irecv (ubuffer , 1 , rowhalo , sourceup , 124 , new , &request[2]) ;
    }

    if (destdown != MPI_PROC_NULL){
        MPI_Send (&local_M[local_dim[0] * local_dim[1]  - kdim * local_dim[0]] , 1 , rowhalo , destdown , 124 , new) ;
    }

    if(sourceup != MPI_PROC_NULL){
        MPI_Wait ( &request[2] , &status[2]) ;
    }

    int sourcedown , destup ;

    MPI_Cart_shift (new , 1 , -1 , &sourcedown , &destup) ;

    if (sourcedown != MPI_PROC_NULL){
        MPI_Irecv (dbuffer , 1 , rowhalo , sourcedown , 125 , new , &request[3]) ;
    }

    if (destup != MPI_PROC_NULL){
        MPI_Send (local_M , 1 , rowhalo , destup , 125 , new) ;
    }

    if(sourcedown != MPI_PROC_NULL){
        MPI_Wait ( &request[3] , &status[3]) ;
    } 
    /*for (int disp = -1 ; disp < 2 ; disp = disp + 2){
        for (int dir = 0 ; dir < 2 ; dir++){
            MPI_Cart_shift (new , dir , disp , &source , &dest) ;

            if (source != MPI_PROC_NULL && disp == -1 && dir == 0){
                MPI_Irecv (rbuffer , 1 , colhalo , source , 123 , new , &request) ; 
            }

            if (source != MPI_PROC_NULL && disp == 1 && dir == 0){
                MPI_Irecv (lbuffer , 1 , colhalo , source , 123 , new , &request) ;
            }

            if (source != MPI_PROC_NULL && disp == -1 && dir == 1){
                MPI_Irecv (dbuffer , 1 , rowhalo , source , 123 , new , &request) ;
            }

            if (source != MPI_PROC_NULL && disp == 1 && dir == 1){
                MPI_Irecv (ubuffer , 1 , rowhalo , source , 123 , new , &request) ;
            }

            if (dest != MPI_PROC_NULL && disp == -1 && dir == 0){
                MPI_Send (local_M , 1 , colhalo , dest , 123 , new ) ;
            }

            if (dest != MPI_PROC_NULL && disp == 1 && dir == 0){
                MPI_Send (&local_M[local_dim[0] - kdim] , 1 , colhalo , dest , 123 , new ) ;
            }

            if (dest != MPI_PROC_NULL && disp == -1 && dir == 1){
                MPI_Send (local_M , 1 , rowhalo , dest , 123 , new ) ;
            }

            if (dest != MPI_PROC_NULL && disp == 1 && dir == 1){
                MPI_Send (&local_M[(local_dim[0] * local_dim[1] - kdim * local_dim[0])] , 1 , rowhalo , dest , 123 , new ) ;
            }
            MPI_Wait( &request , &status) ;
        }
    }*/


    MPI_Barrier (new) ;
   
    MPI_Type_free (&colhalo) ;
    MPI_Type_free (&rowhalo) ;
    MPI_Type_free (&resizedtype) ;
    free (local_M) ;
    free (ubuffer) ;
    free (dbuffer) ;
    free (lbuffer) ;
    free (rbuffer) ;

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