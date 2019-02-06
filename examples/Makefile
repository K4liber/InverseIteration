# AMGX build path
AMG_PATH=/scratch/dteam002# edit this before compile
AMGX_LIB_PATH=${AMG_PATH}/AMGX/build
AMGX_INCLUDE_PATH= ${AMG_PATH}/AMGX/base/include
# InverseIterator build path
II_PATH=/home/dteam002/project/InverseIteration/MPI
# Cuda path
CUDA_PATH=/usr/local/cuda

CFLAGS= -O3 -gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_37,code=sm_37 -gencode=arch=compute_60,code=sm_60 -lm
MPI= -L/opt/openmpi-gcc721-Cuda90/3.1.1/lib -lmpi

hamiltonian:
	mpic++ ${CFLAGS} hamiltonian.cpp -o hamiltonian \
		${MPI} \
		-L${AMGX_LIB_PATH} -lamgxsh \
		-L${II_PATH} -lInverseIterator \
		-L${CUDA_PATH}/lib64 -lcudart \
		-I${II_PATH} \
		-I${CUDA_PATH}/include \
		-I${AMGX_INCLUDE_PATH}

clean:
	rm -f hamiltonian