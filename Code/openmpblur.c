#include <stdio.h>
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



void kernel (double * K , int dim_kernel , int type_kernel , double weight) ;
void write_pgm_image( void *image, int maxval, int xsize, int ysize, const char *image_name);
void read_pgm_image( void **image, int *maxval, int *xsize, int *ysize, const char *image_name);
void swap_image( void *image, int xsize, int ysize, int maxval );
void conv( unsigned short int * result_M , unsigned short int * temp , double * K , int kdim , int * local_dim , int dim_kernel , int i , int j) ;

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

    
    double start_time , end_time ;

    start_time = omp_get_wtime() ;

    read_pgm_image((void **)&M , &maxval , &xsize , &ysize , image_name) ;
    int N[2] = {xsize , ysize} ;
    swap_image((void*)M , xsize , ysize , maxval);
    
    
    
    int type_kernel ;
    int dim_kernel ;
    type_kernel = atoi (argv[2]) ;
    dim_kernel = atoi (argv[3]) ;
    float weight = 0;
    int index = 4;
    if(type_kernel == 1){
    weight = atof (argv[4]) ;
    index = 5;
    }


    int kdim = dim_kernel / 2 ;
    double * K ;
    K = (double*)malloc( dim_kernel * dim_kernel * sizeof(double)) ;
    kernel ( K , dim_kernel , type_kernel , weight) ;

    unsigned short int *temp;
    temp = (unsigned short int*)calloc((N[0] + 2 * kdim)*(N[1] + 2 * kdim), sizeof(unsigned short int));
    unsigned short int * new_M;
    
  

    #pragma omp parallel 
        {
        int i , j;
        #pragma omp single
        {
        int nthreads = omp_get_num_threads() ;
        printf("Number of threads : %d\n" , nthreads) ;
        }
        #pragma omp for schedule(static)
        for (int k = 0 ; k < N[0] * N[1] ; k++){
          i = k / N[0] ;
          j = k % N[0] ;
          temp[(i + kdim) * (N[0]  + 2 * kdim) + j + kdim] = M[i * N[0] + j] ;
        }
        
           #pragma omp for schedule(static)


                for (int k = 0 ; k < N[0] * N[1] ; k++){
                    i = k / N[0] ;
                    j = k % N[0] ;
                    conv(M , temp , K , kdim , N , dim_kernel , i , j) ;
                }

        }
    swap_image((void *)M , N[0] , N[1] , maxval);

    const char * new_image = argv[index] ;
    write_pgm_image((void*)M , maxval , N[0]  , N[1] , new_image) ;

    free(M) ;
    free(temp) ;
    free(K) ;

    end_time = omp_get_wtime() ;
    printf(" Walltime: %10.8f\n" , end_time - start_time) ;
    return 0;
}


void kernel (double * K, int dim_kernel , int type_kernel , double weight)
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
              if( (i == dim_kernel/2) && (j == i))
                K[i*dim_kernel + j] = weight ;
              else  
                    K[i*dim_kernel + j] = (1 - weight)/(float)(dim_kernel*dim_kernel -1) ;
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




void conv( unsigned short int * result_M , unsigned short int * temp , double * K , int kdim , int * local_dim , int dim_kernel , int i , int j )
{
  double pixel = 0;
  int u , v ;
  result_M[i * local_dim[0] + j] = 0;
            #pragma simd
               for(u = 0 ; u < dim_kernel ; u++){
              #pragma simd
                 for(v = 0 ; v < dim_kernel ; v++){
                   
                   //if ( (i + kdim - u) < 0 || (i + kdim - u) >= local_dim[1] || (j + v - kdim) < 0 || (j + v - kdim) >= local_dim[0])
                   // continue;
                   //else
                   pixel += (temp[ (i + kdim + (kdim - u)) * (local_dim[0] + 2 * kdim) + (j + kdim + (v - kdim))] *  K[(dim_kernel - 1 - u) * dim_kernel + v]) ;
                 }
               }
               result_M[i * local_dim[0] + j] = pixel;
            
}