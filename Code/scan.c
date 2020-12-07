
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

int main( int argc , char* argv[]){
    
    int a[6] = {9 , 5 , 1 , 12 , 3 , 7} ;
    int sum[6] = {0 , 0 , 0 , 0 , 0 , 0} ;
    struct  timespec ts;
    double startTime , endTime ;

    startTime = CPU_TIME ;
    #if defined(_OPENMP)

    int nthreads ;
    int s = 0 ;
        #pragma omp single 
        {
        nthreads = omp_get_num_threads() ;
        }

        #pragma omp parallel for 
        for (int i = 0 ; i < 6 ; i ++){
            #pragma omp parallel for reduction (+:sum[i])
            for ( int k = 0 ; k < i + 1 ; k ++ ){
            
            sum[i] += a[k] ;  
            }
              
        }
    
        

        #pragma omp single
        {
            for(int i = 0 ; i < 6 ; i++)
                printf("%d\n" , sum[i]) ;
        endTime = CPU_TIME ;
        printf("Walltime omp : %g seconds\n" , endTime - startTime ) ;
        }

    #else
    sum[0] = a[0];
    printf("%d\n" , sum[0]) ;
    for (int i = 1 ; i < 6 ; i ++){
        sum[i] = a[i] + sum[i -1] ;
        printf("%d\n" , sum[i]) ;
    }

    endTime = CPU_TIME ;
    printf("Walltime : %g seconds\n" , endTime - startTime ) ;

    #endif
    
    return 0 ;
}