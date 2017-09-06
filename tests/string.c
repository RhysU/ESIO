//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// ESIO 0.2.0: ExaScale IO library for turbulence simulation restart files
// http://github.com/RhysU/ESIO
//
// Copyright (C) 2010-2017 The PECOS Development Team
//
// This file is part of ESIO.
//
// ESIO is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3.0 of the License, or
// (at your option) any later version.
//
// ESIO is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with ESIO.  If not, see <http://www.gnu.org/licenses/>.
//
//-----------------------------------------------------------------------el-
// $Id$

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <esio/esio.h>
#include <esio/error.h>

#include "testutils.h"

// Include FCTX and silence useless warnings
#ifdef __INTEL_COMPILER
#pragma warning(push,disable:981)
#endif
#include "fct.h"
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

// Add command line options
static const fctcl_init_t my_cl_options[] = {
    {
        "--preserve",
        "-t",
        FCTCL_STORE_TRUE,
        "Are temporary filenames displayed and files preserved?"
    },
    FCTCL_INIT_NULL /* Sentinel */
};

FCT_BGN()
{
    int preserve = 0;

    // MPI setup: MPI_Init and atexit(MPI_Finalize)
    int world_size, world_rank;
    MPI_Init(&argc, &argv);
    atexit((void (*) ()) MPI_Finalize);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Install the command line options defined above.
    fctcl_install(my_cl_options);

    // Retrieve and sanity check problem size options
    preserve = fctcl_is("--preserve");

    // Obtain default HDF5 error handler
    H5E_auto2_t hdf5_handler;
    void *hdf5_client_data;
    H5Eget_auto2(H5E_DEFAULT, &hdf5_handler, &hdf5_client_data);

    // Obtain default ESIO error handler
    esio_error_handler_t * const esio_handler = esio_set_error_handler_off();
    esio_set_error_handler(esio_handler);

    // Fixture-related details
    const char * const input_dir  = getenv("ESIO_TEST_INPUT_DIR");
    const char * const output_dir = getenv("ESIO_TEST_OUTPUT_DIR");
    char * filetemplate = create_testfiletemplate(output_dir, __FILE__);
    (void) input_dir;  // Possibly unused
    char * filename = NULL;
    esio_handle state;

    FCT_FIXTURE_SUITE_BGN(line)
    {
        FCT_SETUP_BGN()
        {
            ESIO_MPICHKQ(MPI_Barrier(MPI_COMM_WORLD)); // Synchronize

            // Restore HDF5/ESIO default error handling
            H5Eset_auto2(H5E_DEFAULT, hdf5_handler, hdf5_client_data);
            esio_set_error_handler(esio_handler);

            // Rank 0 generates a unique filename and broadcasts it
            int filenamelen;
            if (world_rank == 0) {
                filename = create_testfilename(filetemplate);
                if (preserve) {
                    printf("\nfilename: %s\n", filename);
                }
                filenamelen = strlen(filename);
            }
            ESIO_MPICHKR(MPI_Bcast(&filenamelen, 1, MPI_INT,
                                   0, MPI_COMM_WORLD));
            if (world_rank > 0) {
                filename = calloc(filenamelen + 1, sizeof(char));
            }
            ESIO_MPICHKR(MPI_Bcast(filename, filenamelen, MPI_CHAR,
                                    0, MPI_COMM_WORLD));
            fct_req(filename);

            // Initialize ESIO state
            state = esio_handle_initialize(MPI_COMM_WORLD);
            fct_req(state);
        }
        FCT_SETUP_END();

        FCT_TEARDOWN_BGN()
        {
            // Finalize ESIO state
            esio_handle_finalize(state);

            // Clean up the unique file and filename
            if (world_rank == 0) {
                if (!preserve) unlink(filename);
            }
            free(filename);
        }
        FCT_TEARDOWN_END();

        // Test scalar-valued lines, including overwrite details
        FCT_TEST_BGN(strings)
        {
            const char * msg1 = "Tra la la";
            const char * msg2 = "The quick brown fox jumps over the lazy dog";

            // Open file
            fct_req(0 == esio_file_create(state, filename, 1));

            // Write empty message to disk
            fct_req(0 == esio_string_set(state, "/", "msg1", ""));

            // Overwrite message with test data
            fct_req(0 == esio_string_set(state, "/", "msg1", msg1));

            // Write a second message to disk
            fct_req(0 == esio_string_set(state, "/", "msg2", msg2));

            // Retrieve the first and ensure it comes back cleanly
            {
                char *retrieved1 = esio_string_get(state, "/", "msg1");
                fct_req(retrieved1);
                fct_chk_eq_str(retrieved1, msg1);
                fct_chk(retrieved1 != msg1);
                free(retrieved1);
            }

            // Retrieve the second and ensure it comes back cleanly
            {
                char *retrieved2 = esio_string_get(state, "/", "msg2");
                fct_req(retrieved2);
                fct_chk_eq_str(retrieved2, msg2);
                fct_chk(retrieved2 != msg2);
                free(retrieved2);
            }

            // Attempt to retrieve a nonexistent third and ensure NULL result
            fct_chk(0 == esio_string_get(state, "/", "nonexistent"));

            // Close the file
            fct_req(0 == esio_file_close(state));

            // Reopen the file using normal HDF5 APIs on root processor
            // Ensure we can retrieve the data by other means
            if (world_rank == 0) {
                char *retrieved1 = malloc((strlen(msg1) + 1));
                fct_req(retrieved1);
                const hid_t file_id
                    = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
                fct_req(0 <= H5LTget_attribute_string(file_id, "/",
                                                      "msg1", retrieved1));
                fct_chk_eq_str(retrieved1, msg1);
                free(retrieved1);
                fct_req(0 <= H5Fclose(file_id));
            }
        }
        FCT_TEST_END();
    }
    FCT_FIXTURE_SUITE_END();

    free(filetemplate);
}
FCT_END()
