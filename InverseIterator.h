#ifndef _InverseIterator_H_
#define _InverseIterator_H_

#include "amgx_c.h"

/** 
    Class InverseIterator is intended to perform inverse iteration 
    algorithm to find the smallest eigenvalue of a given matrix,
    using AMGX library (https://github.com/NVIDIA/AMGX) to solve 
    the system of linear equations.
 */
class InverseIterator {

public:

    InverseIterator(double**, int, double, char*);
    double getEigenValue(bool);

private:

    AMGX_matrix_handle A;
    AMGX_vector_handle b, x;
    AMGX_config_handle cfg;
    AMGX_solver_handle solver;
    AMGX_resources_handle res = NULL;
    AMGX_Mode mode = AMGX_mode_dDDI;
    double* h_b;
    double* h_x;
    double** h_matrix;
    double epsilon;
    int N;

    void saveMatrixAsMTX();
    void normalizeSolution();
    double getResiduum();
};

#endif