# AMGX build path
AMG_PATH=/scratch/dteam002 # edit this before compile
AMGX_LIB_PATH=${AMG_PATH}/AMGX/build
AMGX_INCLUDE_PATH= ${AMG_PATH}/AMGX/base/include
CFLAGS= -fPIC -shared -std=c++11

libInverseIterator.so:
	gcc ${CFLAGS} InverseIterator.cpp -o libInverseIterator.so \
		-L${AMGX_LIB_PATH} -lamgxsh -I${AMGX_INCLUDE_PATH}

clean:
	rm -f libInverseIterator.so