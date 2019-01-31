#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "amgx_c.h"

using namespace std;

class InverseIterator {

    public:

        InverseIterator(String configFileName, String matrixFileName, int N) {
            /* init */
            AMGX_SAFE_CALL(AMGX_initialize());
            AMGX_SAFE_CALL(AMGX_initialize_plugins());
            /* system */
            AMGX_SAFE_CALL(AMGX_install_signal_handler());
            //Read config file
            AMGX_config_create_from_file(&cfg, configFileName);

            //Create resources based on config
            AMGX_resources_create_simple(&res, cfg);

            //Create solver object, A,x,b, set precision
            AMGX_solver_create(&solver, res, mode, cfg);
            AMGX_matrix_create(&A,res,mode);
            AMGX_vector_create(&x,res,mode); 
            AMGX_vector_create(&b,res,mode);

            //Read coefficients from a file    
            AMGX_SAFE_CALL(AMGX_read_system(A, x, b, argv[2]));
            int n = 0;
            int xsize_x = 0, xsize_y = 0;
            AMGX_SAFE_CALL(AMGX_matrix_get_size(A, &n, &xsize_x, &xsize_y));
            AMGX_SAFE_CALL(AMGX_vector_set_random(x, n));

            //Try to pass something to the solution vector
            int N = 2;
            *h_b = malloc(N * sizeof(double));
            *h_x = malloc(N * sizeof(double));
            h_b[0] = 5.5;
            h_b[1] = 7;
            AMGX_vector_upload(b, N, 1, h_b);
        }

        double getEigenValue() {

        }

    private:

        AMGX_matrix_handle A;
        AMGX_vector_handle b, x;
        AMGX_config_handle cfg;
        AMGX_solver_handle solver;
        AMGX_resources_handle res = NULL;
        AMGX_Mode mode = AMGX_mode_dDDI;
        double *h_b;
        double *h_x;

        /* *** Save matrix of class Eigen::MatrixXd to mtx format *** */
        void saveMatrixAsMTX(Eigen::MatrixXd tab, int n){
            ofstream file;
            file.open("matrix.mtx");
            file << '%%MatrixMarket matrix coordinate real general' << '\n';
            int nozeros=0;
            for (int i = 0; i<n; ++i){
                for(int j = 0; j<n; ++j){
                    if( tab(i,j) != 0){
                        nozeros++;
                        file << i+1 <<" "<< j+1 <<" "<< tab(i,j) << '\n';
                    }
                }
            }
            file << n << " " << n << " " << nozeros;	
            file.close();	
        }
}