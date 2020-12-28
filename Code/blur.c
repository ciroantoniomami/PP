#include <stdio.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#if ((0x100 & 0xf) == 0x0)
#define I_M_LITTLE_ENDIAN 1
#define swap(mem) (( (mem) & (short int)0xff00) >> 8) +	\
  ( (mem) & (short int)0x00ff << 8)
#else
#define I_M_LITTLE_ENDIAN 0
#define swap(mem) (mem)
#endif

void kernel (double * K , int dim_kernel , int type_kernel , int weight) ;
void write_pgm_image( void *image, int maxval, int xsize, int ysize, const char *image_name);
void read_pgm_image( void **image, int *maxval, int *xsize, int *ysize, const char *image_name);



int main (int argc , char * argv[]){
    
    //Initialize matrix for the image and the kernel
    if (argc <=1){
        fprintf (stderr , " Required at least 3 parameters : PGM Image ,"
         "type of the kernel , dimension of the kernel and weight (if needed) ") ;
        
        exit (-1) ;
    }
    
    const char *image_name = argv[1] ;
    int  xsize ;
    int  ysize ;
    
    int maxval ;
    short int * M ;
    double * new_M ;
    
    

   // M = (short int*)malloc(N[0] * N[1] * sizeof(short int));
    
    
    
   
    read_pgm_image((void **)&M , &maxval , &xsize , &ysize , image_name) ;
    int N[2] = {xsize , ysize} ;
    new_M = (double*)malloc(N[0] * N[1] * sizeof(double));

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
    MPI_Status status[2] ;
    MPI_Request request[2] ;
    //Create a new comm world
    int size;
    MPI_Comm_size(MPI_COMM_WORLD , &size);

    int dims[2] = {0,0};
    MPI_Dims_create(size ,2 , dims);
  
    int periods[2] = {0, 0};
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

    MPI_Datatype type2, resizedtype2 ;
    

    MPI_Type_create_subarray (2, N , local_dim , starts , MPI_ORDER_C , MPI_DOUBLE , &type2) ;
    MPI_Type_create_resized (type2, 0, local_dim[0] * sizeof(double), &resizedtype2) ;
    MPI_Type_commit (&resizedtype2) ;

    int counts[size] ;
    int counts2[size] ;
    int displs[size] ;
    int displs2[size] ;




    if (rank == 0){


        for (int i = 0 ; i < dim_kernel ; i++){
            for (int j = 0 ; j < dim_kernel ; j++)
                printf (" %9.2f ",K[i * dim_kernel + j]) ;
            printf ("\n") ;
        }

    for (int i = 0; i < size; i++) counts[i] = 1 ;
    for (int i = 0; i < size; i++) counts2[i] = local_dim[0] * local_dim[1] ;
    for (int i = 0; i < size; i++) displs[i] = i % dims[0] + i / dims[0] * local_dim[1] * dims[0] ;
    for (int i = 0; i < size; i++) displs2[i] = i % dims[0] * local_dim[0] + i / dims[0] * local_dim[1] * dims[0] * local_dim[0] ;
    }
    MPI_Bcast (counts , size , MPI_INT , 0 , new) ;
    MPI_Bcast (displs , size , MPI_INT , 0 , new) ;

    MPI_Barrier (new) ;
    /* send submatrices to all processes */
    MPI_Scatterv (M , counts , displs , resizedtype , local_M ,
     local_dim[0] * local_dim[1] , MPI_SHORT, 0, new) ;

    MPI_Barrier (new) ;
    
    if( rank == 3){
        const char * new_image = argv[5] ;
        write_pgm_image(local_M , 65535 , local_dim[0] , local_dim[1] , new_image) ;
    }
    // The structure that will be used to exchange data between   
    short int * ubuffer = (short int *) calloc (  local_dim[0] * kdim , sizeof ( short int)) ;    
    short int * dbuffer = (short int *) calloc (  local_dim[0] * kdim , sizeof ( short int)) ;
    short int * lbuffer = (short int *) calloc ( local_dim[1] * kdim , sizeof ( short int)) ; 
    short int * rbuffer = (short int *) calloc ( local_dim[1] * kdim , sizeof ( short int)) ; 
    short int * topleftbuffer = (short int *) calloc ( kdim * kdim , sizeof ( short int)) ; 
    short int * toprightbuffer = (short int *) calloc ( kdim * kdim , sizeof ( short int)) ; 
    short int * downleftbuffer = (short int *) calloc ( kdim * kdim , sizeof ( short int)) ; 
    short int * downrightbuffer = (short int *) calloc ( kdim * kdim , sizeof ( short int)) ; 

    int subsizes_col[2] = {kdim , local_dim[1]} ;
    int subsizes_row[2] = {local_dim[0] , kdim} ;
    int subsizes_angle[2] = {kdim , kdim} ;
    
    MPI_Datatype type_col , colhalo, type_row , rowhalo , type_angle , anglehalo;

    MPI_Type_create_subarray (2, local_dim , subsizes_col , starts , MPI_ORDER_C , MPI_SHORT , &type_col) ;
    MPI_Type_create_resized (type_col, 0, kdim * sizeof(short int), &colhalo) ;
    MPI_Type_commit (&colhalo) ;

    MPI_Type_create_subarray (2, local_dim , subsizes_row , starts , MPI_ORDER_C , MPI_SHORT , &type_row) ;
    //MPI_Type_create_resized (type_row, 0, local_dim[0] * sizeof(short int), &rowhalo) ;
    MPI_Type_commit (&type_row) ;

    MPI_Type_create_subarray (2, local_dim , subsizes_angle , starts , MPI_ORDER_C , MPI_SHORT , &type_angle) ;
    MPI_Type_create_resized (type_angle, 0, kdim * sizeof(short int), &anglehalo) ;
    MPI_Type_commit (&anglehalo) ;

    int left , right ;
    
    int samma ;
    MPI_Cart_shift (new , 1 , 1 , &left , &right) ;
    if (left != MPI_PROC_NULL){
        MPI_Irecv (lbuffer , 1 , colhalo , left , 123 , new , &request[0]) ;
       samma = 1;
    }

    if (right != MPI_PROC_NULL){
        MPI_Send (&local_M[ local_dim[0] - kdim  ] , 1 , colhalo , right , 123 , new) ;
        samma = 2;
    }

    if(left != MPI_PROC_NULL){
        MPI_Wait ( &request[0] , &status[0]) ;
    }
 //   if(samma == 1) {printf("rank %d has recieved on axis x from left\n", rank);
 //   }
 //   if(samma == 2) {printf("rank %d has send on axis x to right\n", rank);
 //   }
   

    if (right != MPI_PROC_NULL){
        MPI_Irecv (rbuffer , 1 , colhalo , right , 124 , new , &request[0]) ;
        samma = 3;
    }

    if (left != MPI_PROC_NULL){
        MPI_Send (local_M , 1 , colhalo , left , 124 , new) ;
        samma = 4;
    }

    if(right != MPI_PROC_NULL){
        MPI_Wait ( &request[0] , &status[0]) ;
    }


//    if(samma == 3) {printf("rank %d has recieved on axis x from right\n", rank);
//    }
//    if(samma == 4) {printf("rank %d has send on axis x to left\n", rank);
//    }

    int up , down ;
    int top_up[2] = { -2 , -2 } ;
    int top_down[2] = { -2 , -2 } ;
    MPI_Cart_shift (new , 0 , 1 , &up , &down) ;

    if (up != MPI_PROC_NULL){
        MPI_Irecv (ubuffer , kdim * local_dim[0] , MPI_SHORT , up , 124 , new , &request[0]) ;

        MPI_Irecv (top_up , 2 , MPI_INT  , up , 12 , new , &request[1]) ;
       
        samma = 5 ;

    }

    int buf[2] = {left , right} ;
    if (down != MPI_PROC_NULL){
        MPI_Send (&local_M[local_dim[0] * local_dim[1]  - kdim * local_dim[0]] , kdim * local_dim[0] , MPI_SHORT , down , 124 , new) ;
        MPI_Send (buf , 2 , MPI_INT  , down , 12 , new) ;

        
        samma = 6 ;    
    }

    if(up != MPI_PROC_NULL){
        MPI_Waitall (2, request , status) ;
    }

 //   if(samma == 5) {printf("rank %d has recieved on axis y from up\n", rank);
 //   }
 //   if(samma == 6) {printf("rank %d has send on axis y to down\n", rank);
 //   }

    if (down != MPI_PROC_NULL){
        MPI_Irecv (dbuffer , kdim * local_dim[0] , MPI_SHORT  , down , 125 , new , &request[0]) ;
        
        MPI_Irecv (top_down , 2 , MPI_INT  , down , 13 , new , &request[1]) ;

        
        samma = 7 ;
    }

    if (up != MPI_PROC_NULL){
        MPI_Send (local_M , kdim * local_dim[0] , MPI_SHORT , up , 125 , new) ;
        
        MPI_Send (buf , 2 , MPI_INT  , up , 13 , new) ;

        samma = 8 ;
    }

    if(down != MPI_PROC_NULL){
        MPI_Waitall (2 , request , status) ;
    } 

//    if(samma == 7) {printf("rank %d has recieved on axis y from down\n", rank);
//    }
//    if(samma == 8) {printf("rank %d has send on axis y to up\n", rank);
//    }
   

 //   printf("I'm rank : %d and my left is : %d , my right : %d\n my up : %d , my down : %d\n", rank ,left ,right, up , down) ;
 //   printf("my up_left : %d , my up_right : %d\n my down_left : %d , my down_right : %d\n", top_up[0] , top_up[1] , top_down[0] , top_down[1] ) ;

    
    if (top_up[0] != MPI_PROC_NULL){
        MPI_Irecv (topleftbuffer , 1 , anglehalo , top_up[0] , 126 , new , &request[0]) ;
        samma = 10 ;
    }

    if (top_down[1] != MPI_PROC_NULL){
        MPI_Send (&local_M[local_dim[0] * local_dim[1] - (kdim -1) * local_dim[0] - kdim] , 1 , anglehalo , top_down[1] , 126 , new) ;
        samma = 11;
    }

    if (top_up[0] != MPI_PROC_NULL){
        MPI_Wait (&request[0] , &status[0]) ;
    }

    if(samma == 10) {printf("rank %d has recieved from topleft\n", rank);
    }
    if(samma == 11) {printf("rank %d has send to downright\n", rank);
    }

//#####################################

    if (top_up[1] != MPI_PROC_NULL){
        MPI_Irecv (toprightbuffer , 1 , anglehalo , top_up[1] , 126 , new , &request[0]) ;
    }

    if (top_down[0] != MPI_PROC_NULL){
        MPI_Send (&local_M[local_dim[0] * local_dim[1] - kdim * local_dim[0]] , 1 , anglehalo , top_down[0] , 126 , new) ;
    }

    if (top_up[1] != MPI_PROC_NULL){
        MPI_Wait (&request[0] , &status[0]) ;
    }

 //#####################################

    if (top_down[0] != MPI_PROC_NULL){
        MPI_Irecv (downleftbuffer , 1 , anglehalo , top_down[0] , 126 , new , &request[0]) ;
    }

    if (top_up[1] != MPI_PROC_NULL){
        MPI_Send (&local_M[local_dim[0] - kdim] , 1 , anglehalo , top_up[1] , 126 , new) ;
    }

    if (top_down[0] != MPI_PROC_NULL){
        MPI_Wait (&request[0] , &status[0]) ;
    }

//#####################################

    if (top_down[1] != MPI_PROC_NULL){
        MPI_Irecv (downrightbuffer , 1 , anglehalo , top_down[1] , 126 , new , &request[0]) ;
    }

    if (top_up[0] != MPI_PROC_NULL){
        MPI_Send (local_M , 1 , anglehalo , top_up[0] , 126 , new) ;
    }

    if (top_down[1] != MPI_PROC_NULL){
        MPI_Wait (&request[0] , &status[0]) ;
    }

  
    double * result_M ;
    result_M = (double *)calloc(local_dim[0] * local_dim[1] , sizeof(double)) ;


    for (int i = kdim ; i < (local_dim[1] - kdim) ; i++){
        for (int j = kdim ; j < (local_dim[0] - kdim) ; j++){
            for (int u = 0 ; u < dim_kernel ; u++){
                for (int v = 0 ; v < dim_kernel ; v++){
                    result_M[i * local_dim[0] + j] += (double)local_M[ (i + (kdim - u)) * local_dim[0] + (j + (v - kdim))] *  K[(dim_kernel - u) * dim_kernel + (v - dim_kernel)] ;
                }
            }
        }
    }

//    for ( int i = 0 ; i < N[0] ; i ++){
//            for (int j = 0 ; j < N[1] ; j++){
//                printf(" %f ", result_M[i * N[0] + j] );
//            }
//            printf("\n");
//        } 


   // MPI_Gatherv(result_M , 1 , resizedtype2 ,
   //             new_M , counts2 , displs2 ,
   //             MPI_DOUBLE , 0 , new ) ;
//
   // if (rank == 0){
   //     const char * new_image = argv[7] ;
   //     write_pgm_image((short int*)new_M , 65535 , N[0] , N[1] , new_image) ;
   //    
   // }

    MPI_Type_free (&colhalo) ;
    MPI_Type_free (&type_row) ;
    MPI_Type_free (&resizedtype) ;
    MPI_Type_free (&anglehalo) ;
    MPI_Type_free (&resizedtype2) ;
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


void write_pgm_image( void *image, int maxval, int xsize, int ysize, const char *image_name)
/*
 * image        : a pointer to the memory region that contains the image
 * maxval       : either 255 or 65536
 * xsize, ysize : x and y dimensions of the image
 * image_name   : the name of the file to be written
 *
 */
{
  FILE* image_file; 
  image_file = fopen(image_name, "w"); 
  
  // Writing header
  // The header's format is as follows, all in ASCII.
  // "whitespace" is either a blank or a TAB or a CF or a LF
  // - The Magic Number (see below the magic numbers)
  // - the image's width
  // - the height
  // - a white space
  // - the image's height
  // - a whitespace
  // - the maximum color value, which must be between 0 and 65535
  //
  // if he maximum color value is in the range [0-255], then
  // a pixel will be expressed by a single byte; if the maximum is
  // larger than 255, then 2 bytes will be needed for each pixel
  //

  int color_depth = 1 + ( maxval > 255 );

  fprintf(image_file, "P5\n# generated by\n# put here your name\n%d %d\n%d\n", xsize, ysize, maxval);
  
  // Writing file
  fwrite( image, 1, xsize*ysize*color_depth, image_file);  

  fclose(image_file); 
  return ;

  /* ---------------------------------------------------------------

     TYPE    MAGIC NUM     EXTENSION   COLOR RANGE
           ASCII  BINARY

     PBM   P1     P4       .pbm        [0-1]
     PGM   P2     P5       .pgm        [0-255]
     PPM   P3     P6       .ppm        [0-2^16[
  
  ------------------------------------------------------------------ */
}


void read_pgm_image( void **image, int *maxval, int *xsize, int *ysize, const char *image_name)
/*
 * image        : a pointer to the pointer that will contain the image
 * maxval       : a pointer to the int that will store the maximum intensity in the image
 * xsize, ysize : pointers to the x and y sizes
 * image_name   : the name of the file to be read
 *
 */
{
  FILE* image_file; 
  image_file = fopen(image_name, "r"); 

  *image = NULL;
  *xsize = *ysize = *maxval = 0;
  
  char    MagicN[2];
  char   *line = NULL;
  size_t  k, n = 0;
  
  // get the Magic Number
  k = fscanf(image_file, "%2s%*c", MagicN );

  // skip all the comments
  k = getline( &line, &n, image_file);
  while ( (k > 0) && (line[0]=='#') )
    k = getline( &line, &n, image_file);

  if (k > 0)
    {
      k = sscanf(line, "%d%*c%d%*c%d%*c", xsize, ysize, maxval);
      if ( k < 3 )
	fscanf(image_file, "%d%*c", maxval);
    }
  else
    {
      *maxval = -1;         // this is the signal that there was an I/O error
			    // while reading the image header
      free( line );
      return;
    }
  free( line );
  
  int color_depth = 1 + ( *maxval > 255 );
  unsigned int size = *xsize * *ysize * color_depth;
  
  if ( (*image = (char*)malloc( size )) == NULL )
    {
      fclose(image_file);
      *maxval = -2;         // this is the signal that memory was insufficient
      *xsize  = 0;
      *ysize  = 0;
      return;
    }
  
  if ( fread( *image, 1, size, image_file) != size )
    {
      free( image );
      image   = NULL;
      *maxval = -3;         // this is the signal that there was an i/o error
      *xsize  = 0;
      *ysize  = 0;
    }  

  fclose(image_file);
  return;
}





