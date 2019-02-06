# AMGX build path
AMG_PATH=/scratch/dteam002# edit this before compile
AMGX_LIB_PATH=${AMG_PATH}/AMGX/build
AMGX_INCLUDE_PATH= ${AMG_PATH}/AMGX/base/include
CFLAGS= -O3 -gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_37,code=sm_37 -gencode=arch=compute_60,code=sm_60 -lm -shared -std=c++11
MPI= -I/opt/openmpi-gcc721-Cuda90/3.1.1/include -L/opt/openmpi-gcc721-Cuda90/3.1.1/lib -lmpi
CUDA_PATH=/usr/local/cuda
X_COMPILER="-fPIC -pthread -Wl,-rpath -Wl,/opt/openmpi-gcc721-Cuda90/3.1.1/lib -Wl,--enable-new-dtags"

libInverseIterator.so:
	nvcc ${CFLAGS} InverseIterator.cpp -o libInverseIterator.so \
		${MPI} \
		-Xcompiler ${X_COMPILER} \
		-L${AMGX_LIB_PATH} -lamgxsh -I${AMGX_INCLUDE_PATH} \
		-L${CUDA_PATH}/lib64 -lcudart \
		-I${CUDA_PATH}/include \

clean:
	rm -f libInverseIterator.so