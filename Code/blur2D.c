#include <stdio.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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

void kernel (double * K , int dim_kernel , int type_kernel , int weight) ;
void write_pgm_image( void *image, int maxval, int xsize, int ysize, const char *image_name);
void read_pgm_image( void **image, int *maxval, int *xsize, int *ysize, const char *image_name);
void swap_image( void *image, int xsize, int ysize, int maxval );
void disorder( unsigned short int * M , unsigned short int * new_M , int * local_dim , int * N , int size , int * dims);
void order( unsigned short int * M , unsigned short int * new_M , int * local_dim , int * N , int size , int * dims);


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

   

    
    if (rank == 3)
    {
        read_pgm_image((void **)&M , &maxval , &xsize , &ysize , image_name) ;
        int N[2] = {xsize , ysize} ;
        int local_dim[2];
    
        local_dim[0] = N[0] / dims[0] ;
        local_dim[1] = N[1] / dims[1] ;
        new_M = (unsigned short int *)malloc(xsize * ysize * sizeof(unsigned short int)) ;
        swap_image((void*)M , xsize , ysize , maxval);


        disorder( M , new_M , local_dim , N , size , dims);
    }

    MPI_Bcast (&xsize , 1 , MPI_INT , 3 , MPI_COMM_WORLD) ;
    MPI_Bcast (&ysize , 1 , MPI_INT , 3 , MPI_COMM_WORLD) ;
    MPI_Bcast (&maxval , 1 , MPI_INT , 3 , MPI_COMM_WORLD) ;

    MPI_Barrier (MPI_COMM_WORLD) ;

    int N[2] = {xsize , ysize} ;
    int local_dim[2];

    local_dim[0] = N[0] / dims[0] ;
    local_dim[1] = N[1] / dims[1] ;

    unsigned short int * local_M ;
    local_M = (unsigned short int *)calloc(local_dim[0] * local_dim[1] , sizeof(unsigned short int)) ;


    int counts[size] ;
    int displs[size] ;

    if (rank == 0){


        for (int i = 0 ; i < dim_kernel ; i++){
            for (int j = 0 ; j < dim_kernel ; j++)
                printf (" %9.2f ",K[i * dim_kernel + j]) ;
            printf ("\n") ;
        }

    for (int i = 0; i < size; i++) counts[i] = local_dim[0] * local_dim[1] ;
    for (int i = 0; i < size; i++) displs[i] = i * local_dim[0] * local_dim[1] ;
    }
    MPI_Bcast (counts , size , MPI_INT , 0 , MPI_COMM_WORLD) ;
    MPI_Bcast (displs , size , MPI_INT , 0 , MPI_COMM_WORLD) ;

    MPI_Barrier (MPI_COMM_WORLD) ;
    /* send submatrices to all processes */
    MPI_Scatterv (new_M , counts , displs , MPI_UNSIGNED_SHORT , local_M ,
     local_dim[0] * local_dim[1] , MPI_UNSIGNED_SHORT, 3, MPI_COMM_WORLD) ;


    //if( rank == 0){
    //    swap_image((void *)local_M , local_dim[0] , local_dim[1] , maxval);
    //
    //    const char * new_image = argv[5] ;
    //    write_pgm_image((void*)local_M , maxval , local_dim[0]  , local_dim[1] , new_image) ;
    //}

    MPI_Gatherv(local_M , local_dim[0] * local_dim[1] , MPI_UNSIGNED_SHORT,
               new_M , counts , displs ,
             MPI_UNSIGNED_SHORT , 3 , MPI_COMM_WORLD ) ;
    if (rank == 3){
        order( M , new_M , local_dim , N , size , dims);
        swap_image((void *)M , N[0] , N[1] , maxval);

        const char * new_image = argv[5] ;
        write_pgm_image((void*)M , maxval , N[0]  , N[1] , new_image) ;
    }
    

    MPI_Finalize() ;
return 0;
}

void disorder( unsigned short int * M , unsigned short int * new_M , int * local_dim , int * N , int size , int * dims){

        for (int k = 0 ; k < size ; k++){
            int u = k / dims[0] ;
            int v = k % dims[0] ;
            for(int i = 0 ; i < local_dim[1] ; i ++){
                for(int j = 0 ; j < local_dim[0] ; j++){
                    new_M[(i + k * local_dim[1]) * local_dim[0] + j] = M[(i + u * local_dim[1]) * N[0] + (j + v * local_dim[0])] ;
                }
            }
        }
 
}

void order( unsigned short int * M , unsigned short int * new_M , int * local_dim , int * N , int size , int * dims){

        for (int k = 0 ; k < size ; k++){
            int u = k / dims[0] ;
            int v = k % dims[0] ;
            for(int i = 0 ; i < local_dim[1] ; i ++){
                for(int j = 0 ; j < local_dim[0] ; j++){
                    M[(i + u * local_dim[1]) * N[0] + (j + v * local_dim[0])] = new_M[(i + k * local_dim[1]) * local_dim[0] + j] ;
                }
            }
        }
} 

void kernel (double * K, int dim_kernel , int type_kernel , int weight)
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
                    K[i*dim_kernel + j] = (1 - weight)/(dim_kernel^2 -1) ;
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

