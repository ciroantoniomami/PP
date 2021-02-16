


all : mpiblur ompiblur openmpblur


mpiblur: mpiblur.c
	mpicc -O3 $< -o $@

ompiblur: mpiblur.c
	mpicc -fopenmp -O3 $< -o $@

openmpblur: openmpblur.c
	gcc -fopenmp -O3 $< -o $@
