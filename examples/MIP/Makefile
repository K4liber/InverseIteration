# AMGX build path
AMG_PATH=/scratch/dteam002# edit this before compile
AMGX_LIB_PATH=${AMG_PATH}/AMGX/build
AMGX_INCLUDE_PATH= ${AMG_PATH}/AMGX/base/include
CFLAGS= -std=c++11

amgx_mpi_capi:
	gcc ${CFLAGS} amgx_mpi_capi.c -o amgx_mpi_capi \
		-L${AMGX_LIB_PATH} -lamgxsh -I${AMGX_INCLUDE_PATH}

clean:
	rm -f amgx_mpi_capi