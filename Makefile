# AMGX build path
AMGX_BUILD_PATH=/home/dteam002/project/AMGX/build
AMGX_INCLUDE_PATH= /home/dteam002/project/AMGX/core/include
CFLAGS= -fPIC -shared -std=c++11

libInverseIterator.so:
	gcc ${CFLAGS} InverseIterator.cpp -o libInverseIterator.so \
		-L${AMGX_BUILD_PATH} -lamgxsh -I${AMGX_INCLUDE_PATH}

clean:
	rm -f libInverseIterator.so