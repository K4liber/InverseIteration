#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include "InverseIterator.h"

/** 
    Initialize the AMGX library.
    It set up computing environment based on configuration file:
    allocate memory, load the matrix and vectors.
    @param matrix The square matrix that we are considering
    @param N Size of the matrix
    @param epsilon Responsible for the accuracy of the calculation
    @param AMGXConfigFilePath Path to the AMGX-type configuration file
 */
InverseIterator::InverseIterator(double** matrix, int N, double epsilon, char* AMGXConfigFilePath) {
    this->epsilon = epsilon;
    this->h_matrix = matrix;
    this->N = N;
    this->h_b = (double*)malloc(N * sizeof(double));
    this->h_x = (double*)malloc(N * sizeof(double));
    AMGX_SAFE_CALL(AMGX_initialize());
    AMGX_SAFE_CALL(AMGX_initialize_plugins());
    AMGX_SAFE_CALL(AMGX_install_signal_handler());
    AMGX_config_create_from_file(&cfg, AMGXConfigFilePath);
    AMGX_resources_create_simple(&res, cfg);
    AMGX_solver_create(&solver, res, mode, cfg);
    AMGX_matrix_create(&A,res,mode);
    AMGX_vector_create(&x,res,mode); 
    AMGX_vector_create(&b,res,mode);
    saveMatrixAsMTX();
    AMGX_SAFE_CALL(AMGX_read_system(A, x, b, "matrix.mtx"));
    remove("matrix.mtx");
    AMGX_SAFE_CALL(AMGX_vector_set_random(x, this->N));
    AMGX_SAFE_CALL(AMGX_vector_set_random(b, this->N));
}

/** Count and return eigenvalue.
    @param log If true - write the progress to the standard output
    @return double
 */
double InverseIterator::getEigenValue(bool log) {
    double res = 1.0;
    int i = 0;
    AMGX_vector_download(b, h_b);
    AMGX_solver_setup(solver, A);
    double* h_xHelp = (double*)malloc(N * sizeof(double));
    while ((res>epsilon) && (i<1000)){
        i++;
        AMGX_solver_solve(solver, b, x);
        AMGX_SAFE_CALL(AMGX_vector_download(x, h_xHelp));
        if (h_xHelp[0] != h_xHelp[0]) {
            break;
        } else {
            h_x = h_xHelp;
        }
        normalizeSolution();
        res = getResiduum();
        if (log) std::cout<<"Itaration: "<<i<<", residual: "<<res<<std::endl;
        AMGX_vector_upload(b, N, 1, h_x);
        AMGX_vector_download(b, h_b);
    }

    double bZero = 0.0;
    for (int i = 0; i < this->N; i++) {
        bZero += h_matrix[0][i] * h_x[i];
    }
    return bZero/h_x[0];
}

void InverseIterator::saveMatrixAsMTX(){
    std::ofstream file;
    file.open("matrix.mtx");
    file << "%%MatrixMarket matrix coordinate real general" << '\n';
    std::vector<std::string> lines;
    std::string line;
    int nozeros=0;
    for (int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            if( h_matrix[i][j] != 0){
                nozeros++;
                line.append(std::to_string(i+1)).append(" ").append(std::to_string(j+1));
                line.append(" ").append(std::to_string(h_matrix[i][j]));
                lines.push_back(line);
                line = "";
            }
        }
    }
    file << N << " " << N << " " << nozeros << "\n";
    for (std::vector<std::string>::iterator it = lines.begin(); it != lines.end(); ++it)
        file << *it << '\n';
    file.close();	
}

void InverseIterator::normalizeSolution() {
    double norm = 0.0;

    for (int i = 0; i < N; i++)
        norm += h_x[i] * h_x[i];
    
    for (int i = 0; i < N; i++)
        h_x[i] = h_x[i]/sqrt(norm);
}

double InverseIterator::getResiduum() {
    double resSub = 0.0;
    double resSum = 0.0;
    double sub = 0.0;
    double sum = 0.0;

    for (int i = 0; i < N; i++) {
        sub = h_x[i] - h_b[i];
        resSub += sub*sub;
        sum = h_x[i] + h_b[i];
        resSum += sum*sum;
    }

    if (resSub > resSum)
        return sqrt(resSum);

    return sqrt(resSub);
}