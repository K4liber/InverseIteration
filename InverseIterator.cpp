#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>

#include "InverseIterator.h"

/* Compile: gcc -fPIC -shared InverseIterator.cpp -o InverseIterator.so */
InverseIterator::InverseIterator() {}
InverseIterator::InverseIterator(double** matrix, double epsilon) {
    this->epsilon = epsilon;
    this->h_matrix = matrix;
    /* init */
    AMGX_SAFE_CALL(AMGX_initialize());
    AMGX_SAFE_CALL(AMGX_initialize_plugins());
    /* system */
    AMGX_SAFE_CALL(AMGX_install_signal_handler());
    //Read config file
    char* configFileName = "FGMRES_AGGREGATION.json";
    AMGX_config_create_from_file(&cfg, configFileName);

    //Create resources based on config
    AMGX_resources_create_simple(&res, cfg);

    //Create solver object, A,x,b, set precision
    AMGX_solver_create(&solver, res, mode, cfg);
    AMGX_matrix_create(&A,res,mode);
    AMGX_vector_create(&x,res,mode); 
    AMGX_vector_create(&b,res,mode);

    saveMatrixAsMTX(matrix, N);
    //Read coefficients from a file 
    AMGX_SAFE_CALL(AMGX_read_system(A, x, b, "matrix.mtx"));
    int xsize_x = 0, xsize_y = 0;
    AMGX_SAFE_CALL(AMGX_matrix_get_size(A, &this->N, &xsize_x, &xsize_y));
    AMGX_SAFE_CALL(AMGX_vector_set_random(x, this->N));
    AMGX_SAFE_CALL(AMGX_vector_set_random(b, this->N));

}

double InverseIterator::getEigenValue() {
    double res = 1.0;
    int i = 0;
    AMGX_vector_download(b, h_b);

    while ((res>epsilon) && (i<1000)){
        i++;
        //Setup and Solve
        AMGX_solver_setup(solver, A);
        AMGX_solver_solve(solver, b, x);
        AMGX_SAFE_CALL(AMGX_vector_download(x, h_x));
        normalize(h_x, this->N);
        res = getNormFromSubstract(h_x, h_b, this->N);
        std::cout<<"Itaration: "<<i<<", residual: "<<res<<std::endl;
        AMGX_vector_upload(x, N, 1, h_b);
    }

    double bZero = 0.0;
    for (int i = 0; i < this->N; i++) {
        bZero += h_matrix[0][i] * h_x[i];
    }
    return bZero/h_x[0];
}

void InverseIterator::saveMatrixAsMTX(double** tab, int n){
    std::ofstream file;
    file.open("matrix.mtx");
    file << "%%MatrixMarket matrix coordinate real general" << '\n';
    std::vector<std::string> lines;
    std::string line;
    int nozeros=0;
    for (int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            if( tab[i][j] != 0){
                nozeros++;
                line.append(sts::string(i+1)).append(" ").append(string(j+1));
                line.append(std::string(tab[i][j]);
                lines.push_back(line);
                line = "";
            }
        }
    }
    file << n << " " << n << " " << nozeros;
    for (std::vector<std::string>::iterator it = lines.begin(); it != lines.end(); ++it)
        file << *it << '\n';
    file.close();	
}

void InverseIterator::normalize(double *v, int n) {
    double norm = 0.0;

    for (int i = 0; i < n; i++)
        norm += v[i] * v[i];
    
    for (int i = 0; i < n; i++)
        v[i] = v[i]/norm;
}

double InverseIterator::getNormFromSubstract(double* v1, double* v2, int n) {
    double norm = 0.0;
    for (int i = 0; i < n; i++) {
        double sub = v1[i] - v2[i];
        norm += sub*sub;
    }
    return norm;
}