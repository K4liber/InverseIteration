#include "InverseIterator.h"
#include <iostream>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <mpi.h>

#define CUDA_CHECK(__err) \
do { \
    if (__err != cudaSuccess) { \
        fprintf(stderr, "Fatal error: %s (at %s:%d)\n", cudaGetErrorString(__err), __FILE__, __LINE__); \
        fprintf(stderr, "*** FAILED - ABORTING\n"); \
        exit(EXIT_FAILURE); \
    } \
} while (0)

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
    double epsilon = atof(argv[3]);
    char AMGXConfigFilePath[] = "FGMRES_AGGREGATION.json";
    float elapsed=0;
    cudaEvent_t start, stop;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));
    double **A = createHamiltonian(N, mu);
    InverseIterator invIter = InverseIterator(A, N, epsilon, AMGXConfigFilePath);

    CUDA_CHECK(cudaEventRecord(start, 0));
    double eigenValue = invIter.getEigenValueMPI(true);
    CUDA_CHECK(cudaEventRecord(stop, 0));

    CUDA_CHECK(cudaEventSynchronize (stop) );
    CUDA_CHECK(cudaEventElapsedTime(&elapsed, start, stop) );
    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));
    
    std::cout<<"Eigenvalue: "<<eigenValue<<std::endl;
    std::cout<<"The elapsed time in gpu was "<<std::to_string(elapsed)<<"ms"<<std::endl;
}