#ifndef _InverseIterator_H_
#define _InverseIterator_H_

#include "amgx_c.h"

/* *** InverseIterator *** 
Class using AMGX library to count eigenvalue of given matrix
using inverse iteration algorithm
*/
class InverseIterator {

public:

    InverseIterator(double**, int, double);
    double getEigenValue();

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