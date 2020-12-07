#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <omp.h>

int value(int x, int y, long  w, long h);

int main(){

    long width =  6000;
    long height = 6000;
    int nthreads = 1;
    double start_time, end_time;
    int* v=(int*)malloc(sizeof(int)*(width*height));
    int i,j;

    FILE *fp;
    if((fp=fopen("mandelbrot.ppm","wt"))!=NULL)
            fprintf(fp,"P3\n%ld %ld 255\n",width,height);
          
    start_time = omp_get_wtime();
    #pragma omp parallel
    {
        #pragma omp master 
        {
            nthreads = omp_get_num_threads();
            printf("Working on %d threads\n",nthreads);
        }

        #pragma omp barrier

        int me = omp_get_thread_num();

        printf("Hi my friend, I'm thread n°%d\n",me);
        
        #pragma omp for schedule(dynamic)

                for(long k=0; k < width*height; k++){
                    i=k/height;
                    j=k%height;
                    v[i*width+j]= value(i,j,width,height);
                }
                    
                
            

       /*#pragma omp master 
        {
        
            for(long i = 0; i < width*height; i++)
                fprintf(fp,"50 %d 150\n",v[i]);
        }*/
            
        
       
        
      
            
        
    }
    fclose(fp);
    free(v);
    end_time = omp_get_wtime();
    printf(" Walltime:%g\n",end_time-start_time);

    return 0;
}


    int value(int x, int y, long w, long h){
        int n_iter=0;
        float complex point = (float)x/w-1.5 + ((float)y/h-0.5)*I;

        float complex z = 0 + 0*I;
            while(abs(z)<2 && n_iter<=50){
                z = z*z + point;
                n_iter++;
            }

            if(n_iter<50) return 100;
            else return 0;
    }