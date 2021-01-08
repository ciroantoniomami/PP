#!/bin/bash                                                   
#PBS -l nodes=1:ppn=24                                       
#PBS -l walltime=01:00:00                                     
                                                              
cd $PBS_O_WORKDIR                                             
module load   openmpi/4.0.3/gnu/9.3.0                         
                                                              
                                       
        
for procs in {1..24}; do
    export OMP_NUM_THREADS=${procs};   
  for i in 1 2 3; do                                      
    ./openmpblur earth-large.pgm 1 11 0.2 result.pgm >>strongOMP11  
  done      
done                                                        
  