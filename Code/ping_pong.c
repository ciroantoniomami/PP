#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#include <unistd.h>

int main(int argc,char* argv[]){
    int rank,size;
    MPI_Status status;
    MPI_Request request;
    MPI_Init(&argc,&argv);
    int PING_PONG_LIMINT=atoi(argv[1]);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int ping_pong_count=0;
    int partner_rank=(rank+1)%2;
    while(ping_pong_count<PING_PONG_LIMINT){
        if(rank==ping_pong_count%2){
            ping_pong_count++;
            MPI_Send(&ping_pong_count,1,MPI_INT,partner_rank,0,MPI_COMM_WORLD);
            printf("%d sent and incremented ping_pong_count %d to %d\n", rank, ping_pong_count,
               partner_rank); 
            sleep(1);    
        }
        else{
            MPI_Recv(&ping_pong_count,1,MPI_INT,partner_rank,0,MPI_COMM_WORLD,&status);
            printf("%d received ping_pong_count %d from %d\n",rank,ping_pong_count,partner_rank);
        }
        
    
    }   
    MPI_Finalize();
}