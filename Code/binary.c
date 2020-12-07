#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>


#if defined(_OPENMP)
#define CPU_TIME (clock_gettime( CLOCK_REALTIME, &ts ), (double)ts.tv_sec + \
		  (double)ts.tv_nsec * 1e-9)

#define CPU_TIME_th (clock_gettime( CLOCK_THREAD_CPUTIME_ID, &myts ), (double)myts.tv_sec +	\
		     (double)myts.tv_nsec * 1e-9)

#else

#define CPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + \
		  (double)ts.tv_nsec * 1e-9)
#endif

int binarysearch( int *data , int start , int end , int Key );

int main ( int argc , char * argv[] ){

    long long int N;
    struct  timespec ts;
    double startTime , endTime ;
    printf("Insert length array\n");
    scanf("%lld" , &N );
    
    int *a = (int*)malloc(sizeof(int) * N) ;
    for(int i = 0 ; i < N ; i++)
        a[i] = i ;
    
    int key ;
    printf("Insert key\n");
    scanf("%d" , &key );

    startTime = CPU_TIME;
    #if defined(_OPENMP)

    int nthreads = omp_get_num_threads() ;
    long long int k = N/nthreads ;
    int control = 0;

    #pragma omp parallel
    
    #pragma omp parallel
    {
    int found = -1 ;
    int id = omp_get_thread_num() ;
    printf("%d" , id) ;
       

        #pragma omp parallel for 
        
        for ( int i = 0 ; i < nthreads ; i ++){
        
        long long int start = i * k ;
        long long int end = start + k - 1 ;

        found = binarysearch ( a , start , end , key ) ;
        if (found != -1){
                printf("found! %d \n" , found ) ;    
                    
                #pragma omp cancel for
            }

        #pragma omp cancellation point for
        }  
        
            
        #pragma omp single 
        {
                found = binarysearch ( a , k*nthreads , N-1 , key ) ;
                if (found != -1)
                {
                    printf("found!\n") ;
                    endTime = CPU_TIME ;

                    printf("Walltime : %g\n" , endTime - startTime);
                }

        }
    
    }           
           
            
    
    
        
       
       
    

    endTime = CPU_TIME ;

    printf("Walltime : %g\n" , endTime - startTime);
    #else
    int found ;
    found = binarysearch (a , 0 , N-1 , key);
    if(found != -1)
            printf("found! %d\n" , found);
    
    endTime = CPU_TIME ;

    printf("Walltime : %g\n" , endTime - startTime);

    #endif
     
    free(a) ;
    return 0 ;

}

int binarysearch( int *data , int start , int end , int Key ){

    int register low = start ;
    int register high = end ;
    int register mid ;
    while ( low <= high ){
        mid = ( low + high ) / 2 ;

        if ( data[mid] < Key )
            low = mid + 1 ;

        else if (data[mid] > Key)
            high = mid - 1 ;
        
        else
            return mid ;
        
        
    }
    return -1 ;
}