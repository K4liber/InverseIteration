# AMGX build path
AMGX_BUILD_PATH=/home/dteam002/project/AMGX/build

CFLAGS= -std=c++11

amgLinearSolver:
	gcc amgLinearSolver.c -o amgLinearSolver \
		-L${AMGX_BUILD_PATH} -lamgxsh

matrix:
	g++ ${CFLAGS} matrix.cpp -o matrix \

clean:
	rm -f amgLinearSolver matrix