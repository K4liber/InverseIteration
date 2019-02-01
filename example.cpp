#include "InverseIterator.h"
#include <iostream>

// Compile: g++ example.cpp -o example -L. -lInverseIterator -lamgxsh
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
}

int main() {
    int N = 10;
    double mu = -5;
    double **A = createHamiltonian(N, mu);
    InverseIterator invIter = InverseIterator(A, 0.0001);
    double eigenValue = invIter.getEigenValue();
    std::cout<<"Eigenvalue: "<<eigenValue<<std::endl;
}