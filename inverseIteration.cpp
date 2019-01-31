#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "amgx_c.h"

using namespace std;

class InverseIteration {

    public:

        InverseIteration(String configFileName) {
            /* init AMGX */
            AMGX_SAFE_CALL(AMGX_initialize());
            AMGX_SAFE_CALL(AMGX_initialize_plugins());
            AMGX_SAFE_CALL(AMGX_install_signal_handler());
            /* Read config file */
            AMGX_config_create_from_file(&cfg, argv[1]);
            /* Create resources based on config */
            AMGX_resources_create_simple(&res, cfg);
        }

    private:

        AMGX_matrix_handle A;
        AMGX_vector_handle b, x;
        AMGX_config_handle cfg;
        AMGX_solver_handle solver;
        AMGX_resources_handle res = NULL;
        AMGX_Mode mode = AMGX_mode_dDDI;


}