#include "InverseIterator.h"
#include <iostream>

// Compile: g++ example.cpp -o example -L/home/dteam002/project/AMGX/build -lamgxsh -L. -lInverseIterator
double** createHamiltonian(int N, double mu) {
    double** A = (double**)malloc(N * sizeof(double*));

    for (int i = 0; i < N; i++) 
        A[i] = (double*)malloc(N * sizeof(double));

    for (int i=0; i < N; i++) {
        A[i][i] = -2.0 - mu;
        A[0][N-1] = 1.0;
        A[N-1][0] = 1.0;	
        if (i > 0) A[i][i-1] = 1.0;
        if (i < N-1) A[i][i+1] = 1.0;
    }
    return A;
}

int main(int argc, char** argv) {
    int N = atoi(argv[1]);
    double mu = atof(argv[2]);
    double **A = createHamiltonian(N, mu);
    InverseIterator invIter = InverseIterator(A, N, 0.01);
    double eigenValue = invIter.getEigenValue();
    std::cout<<"Eigenvalue: "<<eigenValue<<std::endl;
}