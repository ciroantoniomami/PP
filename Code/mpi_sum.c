
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#define SEED 1998
int main(int argc,char *argv[]){
    int id,numprocs;
    double start_time, end_time; 
    MPI_Status status;
    int master=0;
    int tag1=12;
    int tag2=123;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&id);
    long long int N=atoll(argv[1]);
    long long int np=N/numprocs;
    int a[N];
    if(id==0){
    int sum=0;
    int temp=0;
    srand48(SEED);
    for (int i=0;i<N;i++){
    a[i]=rand()%10;
    }
    start_time = MPI_Wtime();
    for(int i=0;i<numprocs-1;i++){
        int j=i*np;
        MPI_Send(&a[j], np ,MPI_INT,(i+1), tag2 ,MPI_COMM_WORLD) ;
    }
    for(int proc=1;proc<numprocs;proc++){
        MPI_Recv(&temp,1,MPI_INT,proc,tag1,MPI_COMM_WORLD,&status) ;
        sum += temp;
    }
    for(int i=(numprocs-1)*np;i<N;i++){
        sum+=a[i];
    }
    printf("sum : %d\n",sum);
    end_time=MPI_Wtime();
    printf ( "\n # walltime : %10.8f \n", end_time - start_time ) ;
    }
    else{
        int temp=0;
        int a[np];
        MPI_Recv(&a,np,MPI_INT,master,tag2,MPI_COMM_WORLD,&status);
        for(int i=0;i<np;i++) temp+=a[i];
        MPI_Send(&temp, 1 ,MPI_INT,master, tag1 ,MPI_COMM_WORLD) ;
    }
    
    start_time = MPI_Wtime();
    int sum;
    int* as = (int*)malloc(sizeof(int) * np);
    int temp=0;
    MPI_Scatter(&a,np,MPI_INT,as,np,MPI_INT,0,MPI_COMM_WORLD);
    for(int i=0;i<np;i++) temp+=as[i];
    MPI_Reduce(&temp,&sum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    free(as);
    if(id==0){
        for(int i=(numprocs)*np;i<N;i++)sum+=a[i];
        end_time=MPI_Wtime();
        printf("sum : %d\n",sum);
        printf ( "\n # walltime : %10.8f \n", end_time - start_time ) ;
    }
    MPI_Finalize() ;
}