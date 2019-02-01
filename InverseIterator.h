#ifndef InverseIterator_h
#define InverseIterator_h

#include <stdlib.h>
#include <string.h>
#include "amgx_c.h"

/* *** InverseIterator *** 
Class using AMGX library to count eigenvalue of given matrix
using inverse iteration algorithm
*/
class InverseIterator {

    public:

        /* *** Class contructor *** */
        InverseIterator(String, double**, int, double);

        /* *** Get matrix eigenvalue using AMGX *** */
        double getEigenValue();

    private:

        AMGX_matrix_handle A;
        AMGX_vector_handle b, x;
        AMGX_config_handle cfg;
        AMGX_solver_handle solver;
        AMGX_resources_handle res = NULL;
        AMGX_Mode mode = AMGX_mode_dDDI;
        double *h_b;
        double *h_x;
        double epsilon;

        /* *** Save matrix to mtx format *** */
        void saveMatrixAsMTX(double**, int);

        /*** Normalize vector v***/
        void normalize(double*, int);

        /*** Substract vectors v1 - v2***/
        double getNormFromSubstract(double*, double*, int);
}

#endif