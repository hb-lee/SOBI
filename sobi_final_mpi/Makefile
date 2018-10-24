CC = mpiicc
PROM = sobi
SOURCE = main_mix.c matrix.c
$(PROM): $(SOURCE)
 #-I/home/pac008/intel/3.2-icc2017-dyn/include -L/home/pac008/intel/3.2-icc2017-dyn/lib -Wall
all:
	$(CC) -Ofast -Wall -fasm-blocks -qopt-prefetch=4 -simd -qopt-dynamic-align -I. -mcmodel=medium -qopenmp -xCORE-AVX2 -ftz -g -I/home/pac008/intel/mkl/include -lmkl_rt -L/home/pac008/intel/mkl/lib/intel64 -o $(PROM) $(SOURCE)

clean:
	rm sobi 
