#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "amgx_c.h"

/*** Compile: gcc -o amgLinearSolver.x amgLinearSolver.c -L/home/dteam002/project/AMGX/build -lamgxsh ***/
/*** Run: ***/
int main(int argc, char** argv) {
    /* init */
    AMGX_SAFE_CALL(AMGX_initialize());
    AMGX_SAFE_CALL(AMGX_initialize_plugins());
    /* system */
    AMGX_SAFE_CALL(AMGX_install_signal_handler());
    AMGX_config_handle cfg;
    AMGX_solver_handle solver;
    AMGX_resources_handle res = NULL;
    AMGX_Mode mode = AMGX_mode_dDDI;
    AMGX_matrix_handle A;
    AMGX_vector_handle b, x;
    //Read config file
    AMGX_config_create_from_file(&cfg, argv[1]);

    //Create resources based on config
    AMGX_resources_create_simple(&res, cfg);

    //Create solver object, A,x,b, set precision
    AMGX_solver_create(&solver, res, mode, cfg);
    AMGX_matrix_create(&A,res,mode);
    AMGX_vector_create(&x,res,mode); 
    AMGX_vector_create(&b,res,mode);
    int N = 2;
    //Read coefficients from a file    
    AMGX_SAFE_CALL(AMGX_read_system(A, x, b, argv[2]));
    int xsize_x = 0, xsize_y = 0;
    //AMGX_SAFE_CALL(AMGX_matrix_get_size(A, &N, &xsize_x, &xsize_y));

    //Try to pass something to the solution vector
    double *h_b = malloc(N * sizeof(double));
    double *h_x = malloc(N * sizeof(double));
    h_b[0] = 5.5;
    h_b[1] = 7;
    AMGX_vector_upload(b, N, 1, h_b);

    //Setup and Solve
    AMGX_solver_setup(solver, A);
    AMGX_solver_solve(solver, b, x);

    //Try to get x vector
    AMGX_SAFE_CALL(AMGX_vector_download(x, h_x));
    printf("Result:\n");
    for (int i=0;i<N;i++) {
        printf("%f\n", h_x[i]);
    }
    return 0;
}
