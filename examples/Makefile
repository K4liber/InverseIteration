# AMGX build path
AMGX_BUILD_PATH=/home/dteam002/project/AMGX/build
# InverseIterator build path
II_PATH=/home/dteam002/project/InverseIteration
# Cuda path
CUDA_PATH=/usr/local/cuda

CFLAGS= -O3 -gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_37,code=sm_37 -gencode=arch=compute_60,code=sm_60 -lm
 

hamiltonian:
	nvcc ${CFLAGS} InverseIterator.cpp -o libInverseIterator.so \
		-L${AMGX_BUILD_PATH} -lamgxsh \
		-L${II_PATH} -lInverseIterator \
		-L${CUDA_PATH}/lib64 -lcudart \
		-I${II_PATH} \
		-I${CUDA_PATH}/include

clean:
	rm -f hamiltonian