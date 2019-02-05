# AMGX build path
AMG_PATH=/scratch/dteam002# edit this before compile
AMGX_LIB_PATH=${AMG_PATH}/AMGX/build
AMGX_INCLUDE_PATH= ${AMG_PATH}/AMGX/base/include

CUDA_PATH=/usr/local/cuda

CFLAGS= -O3 -gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_37,code=sm_37 -gencode=arch=compute_60,code=sm_60 -lm
XCOMPILER= -pthread
XLINKER= -L/opt/openmpi-gcc721-Cuda90/3.1.1/lib -lmpi

amgx_mpi_capi:
	nvcc ${CFLAGS} amgx_mpi_capi.c -o amgx_mpi_capi \
		 ${XLINKER} \
		-L${AMGX_LIB_PATH} -lamgxsh -I${AMGX_INCLUDE_PATH} \
		-L${CUDA_PATH}/lib64 -lcudart \
		-I${CUDA_PATH}/include \

clean:
	rm -f amgx_mpi_capi