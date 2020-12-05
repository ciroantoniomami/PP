#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <omp.h>

int value(int x, int y, long  w, long h);

int main(){
    double start_time, end_time;
    long width =  6000;
    long height = 6000;
    int v;

    FILE *fp;
    start_time = omp_get_wtime();
     if((fp=fopen("mandelbrot.ppm","wt"))!=NULL){
            fprintf(fp,"P3\n%d %d 255\n",width,height);
            for(long i = 0; i < width; i++){
              for(long j=0; j < height; j++){
                v=value(i,j,width,height);
                fprintf(fp,"50 %d 150\n",v);
              }

            }

     }
            fclose(fp);
    end_time = omp_get_wtime();
    printf(" Walltime:%g\n",end_time-start_time);
    return 0;
}


 int value(int x, int y, long w, long h){
        int n_iter=0;
        float complex point = (float)x/w + (float)y/h*I;

        float complex z = 0 + 0*I;
            while(abs(z)<2 && n_iter<=50){
                z = z*z + point;
                n_iter++;
            }

            if(n_iter<50) return 100;
            else return 0;
    }