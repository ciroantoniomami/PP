#include <stdio.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#define XWIDTH 256
#define YWIDTH 256
#define MAXVAL 65535


#if ((0x100 & 0xf) == 0x0)
#define I_M_LITTLE_ENDIAN 1
#define swap(mem) (( (mem) & (short int)0xff00) >> 8) +	\
  ( ((mem) & (short int)0x00ff) << 8)
#else
#define I_M_LITTLE_ENDIAN 0
#define swap(mem) (mem)
#endif

void kernel (double * K , int dim_kernel , int type_kernel , float weight) ;
void write_pgm_image( void *image, int maxval, int xsize, int ysize, const char *image_name);
void read_pgm_image( void **image, int *maxval, int *xsize, int *ysize, const char *image_name);
void swap_image( void *image, int xsize, int ysize, int maxval );
void disorder( unsigned short int * M , unsigned short int * new_M , int * local_row , int * local_col  , int * N , int size , int * dims);
void order( unsigned short int * M , unsigned short int * new_M , int * local_row , int * local_col , int * N , int size , int * dims);
void sendrightbuffer(unsigned short int * local_M , unsigned short int * sendrbuffer , int kdim , int * local_dim);
void sendleftbuffer(unsigned short int * local_M , unsigned short int * sendlbuffer , int kdim , int * local_dim) ;
void sendtoprightbuffer (unsigned short int * local_M , unsigned short int * sendtrbuffer , int kdim , int * local_dim) ;
void sendtopleftbuffer (unsigned short int * local_M , unsigned short int * sendtlbuffer , int kdim , int * local_dim) ;
void senddownrightbuffer (unsigned short int * local_M , unsigned short int * senddrbuffer , int kdim , int * local_dim) ;
void senddownleftbuffer (unsigned short int * local_M , unsigned short int * senddlbuffer , int kdim , int * local_dim) ;
void conv( unsigned short int * result_M , unsigned short int * temp , double * K , int kdim , int * local_dim , int dim_kernel , int i , int j) ;
void local_dimension(int * N , int * local_row , int * local_col , int size , int * dims) ;


int main (int argc , char * argv[])
{
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
    unsigned short int * M ;
    unsigned short int * new_M ;
    
    double start_time , end_time ;


    
   
    
    int type_kernel ;
    int dim_kernel ;
    int index = 4;
    type_kernel = atoi (argv[2]) ;
    dim_kernel = atoi (argv[3]) ;
    float weight = 0 ;
    if(type_kernel == 1){
    float weight ;
    weight = atof (argv[4]) ;
    index = 5;
    }
    int kdim = dim_kernel / 2 ;
    double * K ;
    K = (double*)malloc( dim_kernel * dim_kernel * sizeof(double)) ;
    kernel ( K , dim_kernel , type_kernel , weight) ;

    int mpi_provided_thread_level ;
    MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &mpi_provided_thread_level);

    if ( mpi_provided_thread_level < MPI_THREAD_FUNNELED ) {
        printf("a problem arise when asking for MPI_THREAD_FUNNELED level\n"); MPI_Finalize();
        exit( 1 );
    }
    MPI_Status status[2] ;
    MPI_Request request[2] ;

    start_time = MPI_Wtime() ;
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

    int local_row[size] ;
    int local_col[size] ;
    

    
    if (rank == 0)
    {
        read_pgm_image((void **)&M , &maxval , &xsize , &ysize , image_name) ;
        int N[2] = {xsize , ysize} ;

        local_dimension(N , local_row , local_col , size , dims) ; 
  
        new_M = (unsigned short int *)malloc(xsize * ysize * sizeof(unsigned short int)) ;
        swap_image((void*)M , xsize , ysize , maxval);


        disorder( M , new_M , local_row ,local_col , N , size , dims);
    }

    
    
    MPI_Bcast (&xsize , 1 , MPI_INT , 0 , new) ;
    MPI_Bcast (&ysize , 1 , MPI_INT , 0 , new) ;
    MPI_Bcast (&maxval , 1 , MPI_INT ,0 , new) ;
    MPI_Bcast (local_row , size , MPI_INT , 0 ,new) ;
    MPI_Bcast (local_col , size , MPI_INT , 0 ,new) ;

    //printf("row: %d , col : %d \n" , local_row[rank] , local_col[rank]) ;
    MPI_Barrier (new) ;
    int local_dim[2] = {local_col[rank] , local_row[rank]} ;

    int N[2] = {xsize , ysize} ;
    

    unsigned short int * local_M ;
    local_M = (unsigned short int *)malloc(local_col[rank] * local_row[rank] * sizeof(unsigned short int)) ;


    int counts[size] ;
    int displs[size] ;

    if (rank == 0){
        //for (int i = 0 ; i < dim_kernel ; i++){
        //    for (int j = 0 ; j < dim_kernel ; j++)
        //        printf (" %9.2f ",K[i * dim_kernel + j]) ;
        //    printf ("\n") ;
        //}

        
    displs[0] = 0 ;

    for (int i = 0; i < size; i++) counts[i] = local_col[i] * local_row[i] ;
    for (int i = 1; i < size; i++){
     
    
      displs[i] = displs[i-1] + counts[i-1];

    }
    }
    MPI_Bcast (counts , size , MPI_INT , 0 , new) ;
    MPI_Bcast (displs , size , MPI_INT , 0 , new) ;

    MPI_Barrier (new) ;
    /* send submatrices to all processes */
    MPI_Scatterv (new_M , counts , displs , MPI_UNSIGNED_SHORT , local_M ,
     local_col[rank] * local_row[rank] , MPI_UNSIGNED_SHORT, 0, new) ;


    //if( rank == 0){
    //    swap_image((void *)local_M , local_dim[0] , local_row[rank] , maxval);
    //
    //    const char * new_image = argv[5] ;
    //    write_pgm_image((void*)local_M , maxval , local_dim[0]  , local_row[rank] , new_image) ;
    //}

    // The structure that will be used to exchange data between   
    unsigned short int * ubuffer = (unsigned short int *) calloc (  local_col[rank]  * kdim , sizeof (unsigned short int)) ;    
    unsigned short int * downrightbuffer = (unsigned short int *) calloc ( kdim * kdim , sizeof (unsigned short int)) ; 
    unsigned short int * dbuffer = (unsigned short int *) calloc (  local_col[rank]  * kdim , sizeof (unsigned short int)) ;
    unsigned short int * lbuffer = (unsigned short int *) calloc ( local_row[rank] * kdim , sizeof (unsigned short int)) ; 
    unsigned short int * rbuffer = (unsigned short int *) calloc ( local_row[rank] * kdim , sizeof (unsigned short int)) ; 
    unsigned short int * topleftbuffer = (unsigned short int *) calloc ( kdim * kdim , sizeof (unsigned short int)) ; 
    unsigned short int * toprightbuffer = (unsigned short int *) calloc ( kdim * kdim , sizeof (unsigned short int)) ; 
    unsigned short int * downleftbuffer = (unsigned short int *) calloc ( kdim * kdim , sizeof (unsigned short int)) ; 
    unsigned short int * sendlbuffer = (unsigned short int *) calloc ( local_row[rank] * kdim , sizeof (unsigned short int)) ;
    unsigned short int * sendrbuffer = (unsigned short int *) calloc ( local_row[rank] * kdim , sizeof (unsigned short int)) ;
    unsigned short int * sendtrbuffer = (unsigned short int *) calloc ( kdim * kdim , sizeof (unsigned short int)) ; 
    unsigned short int * sendtlbuffer = (unsigned short int *) calloc ( kdim * kdim , sizeof (unsigned short int)) ; 
    unsigned short int * senddrbuffer = (unsigned short int *) calloc ( kdim * kdim , sizeof (unsigned short int)) ; 
    unsigned short int * senddlbuffer = (unsigned short int *) calloc ( kdim * kdim , sizeof (unsigned short int)) ; 


    sendleftbuffer(local_M , sendlbuffer , kdim , local_dim);
    sendrightbuffer(local_M , sendrbuffer , kdim , local_dim);
    sendtoprightbuffer(local_M , sendtrbuffer , kdim ,local_dim) ;
    sendtopleftbuffer(local_M , sendtlbuffer ,kdim , local_dim) ;
    senddownrightbuffer(local_M , senddrbuffer , kdim ,local_dim) ;
    senddownleftbuffer(local_M , senddlbuffer , kdim , local_dim) ;
    

    int starts [2] = {0 , 0} ;

 
    int left , right ;
    
    int samma ;
    MPI_Cart_shift (new , 1 , 1 , &left , &right) ;
  
    MPI_Sendrecv(sendrbuffer , kdim * local_row[rank] , MPI_UNSIGNED_SHORT , right , 123 , 
    lbuffer , kdim * local_row[rank] , MPI_UNSIGNED_SHORT, left , 123 , new , &status[0]) ;
 
   
//#####################################


    MPI_Sendrecv(sendlbuffer , kdim * local_row[rank]  , MPI_UNSIGNED_SHORT , left , 124 , 
    rbuffer , kdim * local_row[rank]  , MPI_UNSIGNED_SHORT , right , 124 , new , &status[0]) ;

//#####################################

    int up , down ;
    int top_up[2] = { MPI_PROC_NULL , MPI_PROC_NULL } ;
    int top_down[2] = { MPI_PROC_NULL , MPI_PROC_NULL } ;
    MPI_Cart_shift (new , 0 , 1 , &up , &down) ;
    int buf[2] = {left , right} ;


    MPI_Sendrecv(&local_M[local_col[rank] * local_row[rank]   - kdim * local_col[rank]] , kdim * local_col[rank]  , MPI_UNSIGNED_SHORT  , down , 124 ,
    ubuffer , kdim * local_col[rank] , MPI_UNSIGNED_SHORT , up , 124 , new , &status[0]) ;

    MPI_Sendrecv(buf , 2 , MPI_INT  , down , 12 ,
    top_up , 2 , MPI_INT  , up , 12 , new , &status[0]) ;


//#####################################  

    MPI_Sendrecv(local_M , kdim * local_col[rank] , MPI_UNSIGNED_SHORT , up , 125 ,
    dbuffer , kdim * local_col[rank] , MPI_UNSIGNED_SHORT  , down , 125 , new , &status[0]) ;
    
    MPI_Sendrecv(buf , 2 , MPI_INT  , up , 13 ,
    top_down , 2 , MPI_INT  , down , 13 , new , &status[0]) ;

//#####################################
   

 //   printf("I'm rank : %d and my left is : %d , my right : %d\n my up : %d , my down : %d\n", rank ,left ,right, up , down) ;
 //   printf("my up_left : %d , my up_right : %d\n my down_left : %d , my down_right : %d\n", top_up[0] , top_up[1] , top_down[0] , top_down[1] ) ;

    
 
    MPI_Sendrecv(senddrbuffer , kdim*kdim , MPI_UNSIGNED_SHORT , top_down[1] , 127 ,
    topleftbuffer , kdim*kdim  , MPI_UNSIGNED_SHORT, top_up[0] , 127 , new , &status[0]) ;

//#####################################


    MPI_Sendrecv(senddlbuffer , kdim*kdim  , MPI_UNSIGNED_SHORT , top_down[0] , 126 , 
    toprightbuffer , kdim*kdim  , MPI_UNSIGNED_SHORT , top_up[1] , 126 , new , &status[0]) ;
   

 //#####################################

  

    MPI_Sendrecv(sendtrbuffer , kdim*kdim , MPI_UNSIGNED_SHORT , top_up[1] , 126 ,
    downleftbuffer , kdim*kdim , MPI_UNSIGNED_SHORT , top_down[0] , 126 , new , &status[0]) ;

//#####################################


    MPI_Sendrecv(sendtlbuffer, kdim*kdim ,MPI_UNSIGNED_SHORT , top_up[0] , 126 ,
    downrightbuffer , kdim*kdim  , MPI_UNSIGNED_SHORT , top_down[1] , 126 , new , &status[0]) ;

    MPI_Barrier(new);
    unsigned short int * temp;
    temp = (unsigned short int*)malloc((local_col[rank] + (2 * kdim)) * (local_row[rank]+ (2 * kdim)) * sizeof(unsigned short int)) ;

    


    #pragma omp parallel 
    {
   

        for (int k = 0 ; k < kdim * local_col[rank]  ; k++){
          int i = k / local_col[rank]  ;
          int j = k % local_col[rank]  ;
          temp[i * (local_col[rank]  + 2 * kdim)  + j + kdim] = ubuffer[i * local_col[rank]  + j] ;
        
        }
      
   
        #pragma omp for schedule(static)
        for (int k = 0 ; k < local_row[rank] * local_col[rank]  ; k++){
          int i = k / local_col[rank] ;
          int j = k % local_col[rank] ;
          temp[(i + kdim) * (local_col[rank]  + 2 * kdim) + j + kdim] = local_M[i * local_col[rank] + j] ;
        }

      
    


      

        for (int k = 0 ; k <  kdim * local_col[rank] ; k++){
          int i = k / local_col[rank] ;
          int j = k % local_col[rank] ;
          temp[(i + local_row[rank] + kdim) * (local_col[rank] + 2 * kdim)  + j + kdim] = dbuffer[i * local_col[rank] + j] ;
        }
     

        for (int k = 0 ; k <  kdim * local_row[rank] ; k++){
          int i = k / kdim ;
          int j = k % kdim ;
          temp[(i + kdim) * (local_col[rank] + 2 * kdim)  + j] = lbuffer[i * kdim + j] ;
        }
      
     

        for (int k = 0 ; k <  kdim * local_row[rank] ; k++){
          int i = k / kdim ;
          int j = k % kdim ;
          temp[(i + kdim) * (local_col[rank] + 2 * kdim)  + j + kdim + local_col[rank]] = rbuffer[i * kdim + j] ;
        }
      
    

        for (int k = 0 ; k <  kdim * kdim ; k++){
          int i = k / kdim ;
          int j = k % kdim ;
          temp[i * (local_col[rank] + 2 * kdim)  + j + kdim + local_col[rank]] = toprightbuffer[i * kdim + j] ;
        }
      
     

        for (int k = 0 ; k <  kdim * kdim ; k++){
          int i = k / kdim ;
          int j = k % kdim ;
          temp[i * (local_col[rank] + 2 * kdim)  + j ] = topleftbuffer[i * kdim + j] ;
        }
      
  

        for (int k = 0 ; k <  kdim * kdim ; k++){
          int i = k / kdim ;
          int j = k % kdim ;
          temp[(i + kdim + local_row[rank]) * (local_col[rank] + 2 * kdim)  + j + kdim + local_col[rank]] = downrightbuffer[i * kdim + j] ;
        }
      
    

        for (int k = 0 ; k <  kdim * kdim ; k++){
          int i = k / kdim ;
          int j = k % kdim ;
          temp[(i + kdim + local_row[rank]) * (local_col[rank] + 2 * kdim)  + j ] = downleftbuffer[i * kdim + j] ;
        }
      
    }

    
      
      
   //if( rank == 1){
   //   swap_image((void *)temp, local_col[rank] + 2 * kdim , local_row[rank] + 2 * kdim , maxval);
   //   const char * new_image = argv[5] ;
   //   write_pgm_image((void*)temp , maxval , local_col[rank]  + 2 * kdim  , local_row[rank] + 2 * kdim, new_image) ;
   //}
   //if( rank == 6){
   //   swap_image((void *)topleftbuffer, kdim , kdim , maxval);
   //   const char * new_image = argv[5] ;
   //   write_pgm_image((void*)topleftbuffer , maxval , kdim , kdim, new_image) ;
   //}


    MPI_Barrier(new);

    unsigned short int * result_M ;
    result_M = (unsigned short int *)malloc(local_col[rank]  * local_row[rank] * sizeof(unsigned short int)) ;
    
  #pragma omp parallel 
  {

    int i , j ;
   #pragma omp for schedule(static)


        for (int k = 0 ; k < local_row[rank] * local_col[rank]  ; k++){
            i = k / local_col[rank] ;
            j = k % local_col[rank] ;
            conv(result_M , temp , K , kdim , local_dim , dim_kernel , i , j) ;
        }
  }
    MPI_Barrier(new);
    MPI_Gatherv(result_M , local_col[rank] * local_row[rank] , MPI_UNSIGNED_SHORT,
               new_M , counts , displs ,
             MPI_UNSIGNED_SHORT , 0 , new ) ;
    
    MPI_Barrier(new);

    if (rank == 0){
        order( M , new_M , local_row ,local_col , N , size , dims);
        swap_image((void *)M , N[0] , N[1] , maxval);

        const char * new_image = argv[index] ;
        write_pgm_image((void*)M , maxval , N[0]  , N[1] , new_image) ;
        free (M) ;
        free (new_M) ;
    }



    free(local_M) ;
    free(ubuffer) ;
    free(dbuffer) ;
    free(lbuffer) ;
    free(rbuffer) ;
    free(toprightbuffer) ;
    free(downrightbuffer) ;
    free(topleftbuffer) ;
    free(downleftbuffer) ;
    free(temp) ;
    free(sendlbuffer) ;
    free(sendrbuffer) ;
    free(sendtrbuffer) ;
    free(sendtlbuffer) ;
    free(senddrbuffer) ;
    free(senddlbuffer) ;





    end_time = MPI_Wtime() ;

    if (rank == 0){
      printf("Executing on %d processors\n" , size) ;
      printf("walltime : %10.8f\n" , end_time - start_time) ;
    }
    MPI_Finalize() ;
    
    free (K) ;

return 0;
}

void disorder( unsigned short int * M , unsigned short int * new_M , int * local_row , int * local_col , int * N , int size , int * dims){

      int sumcol = 0 ;
      int sumrow = 0 ;
      int t = 0 ;
        for (int k = 0 ; k < size ; k++){
            int u = k / dims[1] ;
            int v = k % dims[1] ;
            if (k > 0){
              sumcol += local_col[k - 1] ;
              
              t += local_row[k-1]*local_col[k-1] ;
            }
            if( v == 0 && k > 0){
              sumcol = 0 ;
              sumrow += local_row[k - 1];
            }
            
            for(int i = 0 ; i < local_row[k] ; i ++){
                for(int j = 0 ; j < local_col[k] ; j++){
                    new_M[(t + i * local_col[k]) + j] = M[(i + sumrow) * N[0] + (j +  sumcol)] ;
                }
            }
        }
 
}

void order( unsigned short int * M , unsigned short int * new_M , int * local_row , int * local_col, int * N , int size , int * dims){


        int sumcol = 0 ;
        int sumrow = 0 ;
        
        int t = 0 ;
        for (int k = 0 ; k < size ; k++){
            int u = k / dims[1] ;
            int v = k % dims[1] ;
            if (k > 0){
              sumcol += local_col[k - 1] ;
             
              t += local_row[k-1]*local_col[k-1] ;

            }
            if( v == 0 && k > 0){
              sumcol = 0 ;
              sumrow += local_row[k - 1];
            }
            for(int i = 0 ; i < local_row[k] ; i ++){
                for(int j = 0 ; j < local_col[k] ; j++){
                    M[(i + sumrow) * N[0] + (j + sumcol)] = new_M[(t + i * local_col[k] ) + j] ;
                }
            }
        }
} 

void kernel (double * K, int dim_kernel , int type_kernel , float weight)
{
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
                    K[i*dim_kernel + j] = (1 - weight)/(dim_kernel*dim_kernel -1) ;
            }

        break;
                
    case 2 :

        for(int i = 0 ; i < dim_kernel ; i++)
            for(int j = 0 ; j < dim_kernel ; j++)
                K[i*dim_kernel + j] = 0.2 ;
    
    default :
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

  fprintf(image_file, "P5\n# generated by\n# Ciro Antonio Mami\n%d %d\n%d\n", xsize, ysize, maxval);
  
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

void swap_image( void *image, int xsize, int ysize, int maxval )
/*
 * This routine swaps the endianism of the memory area pointed
 * to by ptr, by blocks of 2 bytes
 *
 */
{
  if ( maxval > 255 )
    {
      // pgm files has the short int written in
      // big endian;
      // here we swap the content of the image from
      // one to another
      //
      unsigned int size = xsize * ysize;
      for ( int i = 0; i < size; i+= 1 )
  	((unsigned short int*)image)[i] = swap(((unsigned short int*)image)[i]);
    }
  return;
}

void sendrightbuffer(unsigned short int * local_M , unsigned short int * sendrbuffer , int kdim , int * local_dim){
    for (int i = 0 ; i < local_dim[1] ; i++){
        for (int j = 0 ; j < kdim ; j++){
            sendrbuffer[i * kdim + j] = local_M[i * local_dim[0] + j + local_dim[0] - kdim] ;
        }
    }
}

void sendleftbuffer(unsigned short int * local_M , unsigned short int * sendlbuffer , int kdim , int * local_dim){
    for (int i = 0 ; i < local_dim[1] ; i++){
        for (int j = 0 ; j < kdim ; j++){
            sendlbuffer[i * kdim + j] = local_M[i * local_dim[0] + j] ;
        }
    }
}

void sendtoprightbuffer (unsigned short int * local_M , unsigned short int * sendtrbuffer , int kdim , int * local_dim){
    for (int i = 0 ; i < kdim ; i++){
        for (int j = 0 ; j < kdim ; j++){
            sendtrbuffer[i * kdim + j] = local_M[i * local_dim[0] + j + local_dim[0] - kdim] ;
        }
    }
}

void sendtopleftbuffer (unsigned short int * local_M , unsigned short int * sendtlbuffer , int kdim , int * local_dim){
    for (int i = 0 ; i < kdim ; i++){
        for (int j = 0 ; j < kdim ; j++){
            sendtlbuffer[i * kdim + j] = local_M[i * local_dim[0] + j] ;
        }
    }
}

void senddownrightbuffer (unsigned short int * local_M , unsigned short int * senddrbuffer , int kdim , int * local_dim){
    for (int i = 0 ; i < kdim ; i++){
        for (int j = 0 ; j < kdim ; j++){
            senddrbuffer[i * kdim + j] = local_M[(i + local_dim[1]- kdim) * local_dim[0] + j + local_dim[0] - kdim] ;
        }
    }
}

void senddownleftbuffer (unsigned short int * local_M , unsigned short int * senddlbuffer , int kdim , int * local_dim){
    for (int i = 0 ; i < kdim ; i++){
        for (int j = 0 ; j < kdim ; j++){
            senddlbuffer[i * kdim + j] = local_M[(i + local_dim[1] - kdim) * local_dim[0] + j] ;
        }
    }
}

void conv( unsigned short int * result_M , unsigned short int * temp , double * K , int kdim , int * local_dim , int dim_kernel , int i , int j )
{
  double pixel = 0;
  int u , v ;
  result_M[i * local_dim[0] + j] = 0;
 // #pragma omp parallel for schedule(static)

                //for (int t = 0 ; t < dim_kernel * dim_kernel ; t++){
                //    u = t / dim_kernel ;
                //    v = t % dim_kernel ;
                //    result_M[i * local_dim[0] + j] += (temp[ ((i + kdim) + (kdim - u)) * (local_dim[0] + 2 * kdim) + (j + kdim + (v - kdim))] *  K[(dim_kernel - 1 - u) * dim_kernel + v]) ;
                //}

                for(u = 0 ; u < dim_kernel ; u++){
                  for(v = 0 ; v < dim_kernel ; v++){
                     pixel += (temp[ ((i + kdim) + (kdim - u)) * (local_dim[0] + 2 * kdim) + (j + kdim + (v - kdim))] *  K[(dim_kernel - 1 - u) * dim_kernel + v]) ;
                  }
                }
                result_M[i * local_dim[0] + j] = pixel;
            
}



void local_dimension(int * N , int * local_row , int * local_col , int size , int * dims){
    
        for(int i = 0 ; i < size ; i++){
            local_row[i] = N[1] / dims[0] ;
            local_col[i] = N[0] / dims[1] ;
        }
        for(int i = (dims[1] - 1) ; i < size ; i+= dims[1]){
            local_col[i] = N[0] - (dims[1] - 1) * ( N[0] / dims[1] ); 

        }

        for(int i = dims[1]*dims[0] - dims[1]; i < size ; i++){
            local_row[i] = N[1] - (dims[0] - 1) * (N[1] / dims[0] ); 

        }
         

    
}