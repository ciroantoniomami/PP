#!/bin/bash                                                   
#PBS -l nodes=1:ppn=48                                       
#PBS -l walltime=01:00:00                                     
                                                              
cd $PBS_O_WORKDIR                                             
module load   openmpi/4.0.3/gnu/9.3.0                         
                                                              
                                       
        
for procs in {1..48}; do
  for i in 1 2 3; do                                         
    mpirun  --mca btl '^openib' -np ${procs} ./mpiblur earth-large.pgm 1 11 0.2 result.pgm >>strongMPI11  
  done      
done                                                        
  
