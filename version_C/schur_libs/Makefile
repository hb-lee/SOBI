CC = gcc
PROM = sobi
SOURCE = main.c light_matrix.c svdcmp.c schur_lib/rt_nonfinite.c schur_lib/rtGetInf.c schur_lib/rtGetNaN.c schur_lib/schur_decompose.c schur_lib/schur_decompose_initialize.c schur_lib/schur_decompose_rtwutil.c schur_lib/schur_decompose_terminate.c schur_lib/xdlanv2.c schur_lib/xgehrd.c schur_lib/xgerc.c schur_lib/xhseqr.c schur_lib/xnrm2.c schur_lib/xrot.c schur_lib/xungorghr.c schur_lib/xzlarfg.c  
$(PROM): $(SOURCE)
	$(CC) -Wall -g -lm -I/home/lhb/PAC/software/mkl/include -lmkl_rt -L/home/lhb/PAC/software/mkl/lib/intel64 -L/home/lhb/PAC/software/lib/intel64 -o $(PROM) $(SOURCE)

