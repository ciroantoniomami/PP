#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#define SEED 91298

int main(int argc,char* argv[])
{

  double start_time, end_time;
  long long unsigned int M = 0;

  int nthreads;
  double pi;
  
  start_time = omp_get_wtime();
  
 #pragma omp parallel  
 #pragma omp master
  {
    nthreads = omp_get_num_threads();    
  }

  long long unsigned int N = atoll(argv[1]);
  printf("omp calculation with %d threads\nN=%Ld\n", nthreads ,N);
  
 #pragma omp parallel
  {
    int myid = omp_get_thread_num();
    int unsigned short myseeds[3] = {SEED+(myid),SEED+(myid*3+1), SEED+(myid*4+2)};
    
    seed48( myseeds );
    
   #pragma omp for reduction(+:M)
    for( long long unsigned int i=0;i<N;i++)
      {
	double x = erand48( myseeds ); 
	double y = erand48( myseeds );
	
	M += ( (x*x + y*y)<=1 );
      } 
  }    
    
  pi = (4.0*(double)M)/N;
  end_time = omp_get_wtime();
  
  printf("Estimation of pi: %1.9f\n Walltime:%g\n",pi,end_time-start_time);
  return 0;
}
