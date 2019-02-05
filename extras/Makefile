# AMGX build path
AMG_PATH=/scratch/dteam002 # edit this before compile
AMGX_LIB_PATH=${AMG_PATH}/AMGX/build
AMGX_INCLUDE_PATH= ${AMG_PATH}/AMGX/base/include
CFLAGS= -std=c++11

amgLinearSolver:
	gcc amgLinearSolver.c -o amgLinearSolver \
		-L${AMGX_LIB_PATH} -lamgxsh -I${AMGX_INCLUDE_PATH}

matrix:
	g++ ${CFLAGS} matrix.cpp -o matrix \

clean:
	rm -f amgLinearSolver matrix