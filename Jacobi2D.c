#include<mpi.h>
#include<stdio.h>
#include <stdlib.h>
#include <stdbool.h>


int main(int argc,char* argv[]){

    MPI_Init(&argc,&argv);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    int dims[2]={0,0};
    MPI_Dims_create(size,2,dims);

    int periods[2]={true,true};
    int reorder=true;

    MPI_Comm new;
    MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,reorder,&new);

    int rank;
    MPI_Comm_rank(new,&rank);

    if(rank==0){
        int N[2];
        printf("Insert problem size on first axis:\n");
        scanf("%d",&N[0]);
        printf("Insert problem size on second axis:\n");
        scanf("%d",&N[1]);
    }

    MPI_Barrier(new);
    
    int up,down,left,right;
    MPI_Cart_shift(new,1,1,&up,&down);
    MPI_Cart_shift(new,0,1,&left,&right);

    int my_coord[2];
    MPI_Cart_coords(new,rank,2,my_coord);
    printf("[Rank: %d] coordinates (%d,%d) up: %d down: %d left: %d right: %d\n",rank,my_coord[0],my_coord[1],up,down,left,right);





    MPI_Finalize();
    return 0;
}