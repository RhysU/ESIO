//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// esio 0.0.1: ExaScale IO library for turbulence simulation restart files
// http://pecos.ices.utexas.edu/
//
// Copyright (C) 2010 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
// $Id$

#if    !defined(TYPE) \
    || !defined(AFFIX)
#error "attribute_template.c should be #included after appropriate #defines"
#endif

// See 'feature_test_macros(7) for details'
#define _GNU_SOURCE

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <esio/esio.h>
#include <esio/error.h>

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
        "--ncomponents",
        "-n",
        FCTCL_STORE_VALUE,
        "Sets the number of components used in vline test"
    },
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
    const int ncomponents = (int) strtol(
        fctcl_val2("--ncomponents","2"), (char **) NULL, 10);
    if (ncomponents < 1) {
        fprintf(stderr, "\n--ncomponents=%d < 1\n", ncomponents);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

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
    (void) input_dir;  // Possibly unused
    (void) output_dir; // Possibly unused
    char * filename = NULL;
    esio_handle state;

    FCT_FIXTURE_SUITE_BGN(attribute_suite)
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
                filename = tempnam(output_dir, "l00t");
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
            assert(filename);

            // Initialize ESIO state
            state = esio_handle_initialize(MPI_COMM_WORLD);
            assert(state);
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
            if (filename) free(filename);
        }
        FCT_TEARDOWN_END();

        // Test scalar-valued attributes, including overwrite details
        FCT_TEST_BGN(attribute)
        {
            TYPE value;

            // Open file
            fct_req(0 == esio_file_create(state, filename, 1));

            // Write zero to disk and flush the buffers
            value = 0;
            fct_req(0 == AFFIX(esio_attribute_write)(state, "attribute", &value));
            fct_req(0 == esio_file_flush(state));

            // Populate value with non-zero data
            // TODO Vary this based on the rank?
            value = 5678;

            // Overwrite zeros on disk with test data
            fct_req(0 == AFFIX(esio_attribute_write)(state, "attribute", &value));

            // Clear storage in memory
            value = 0;

            { // Ensure we can retrieve the size correctly
                int count;
                fct_req(0 == esio_attribute_sizev(state, "attribute", &count));
                fct_chk_eq_int(count, 1);
            }

            // Close the file
            fct_req(0 == esio_file_close(state));

            // Reopen the file using normal HDF5 APIs on root processor
            // Ensure we can retrieve the data by other means
            if (world_rank == 0) {
                TYPE data;
                const hid_t file_id
                    = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
                fct_req(0 <= AFFIX(H5LTget_attribute)(file_id, "/",
                                                      "attribute",
                                                      &data));

                fct_chk(data == 5678);

                fct_req(0 <= H5Fclose(file_id));
            }

            // Re-read the file in a distributed manner and verify contents
            fct_req(0 == esio_file_open(state, filename, 0));
            fct_req(0 == AFFIX(esio_attribute_read)(state, "attribute", &value));
            fct_chk(value == 5678);
            fct_req(0 == esio_file_close(state));
        }
        FCT_TEST_END();

        // Test vector-like attributes
        FCT_TEST_BGN(attribute)
        {
            TYPE *value;

            // Allocate storage
            value = calloc(ncomponents, sizeof(TYPE));
            fct_req(value);

            // Open file
            fct_req(0 == esio_file_create(state, filename, 1));

            // Populate test data
            for (int i = 0; i < ncomponents; ++i) {
                value[i] = (TYPE) i + 5678;
            }

            // Write attribute to file
            fct_req(0 == AFFIX(esio_attribute_writev)(state, "attribute",
                                                      value, ncomponents));

            { // Ensure we can retrieve the size correctly
                int count;
                fct_req(0 == esio_attribute_sizev(state, "attribute", &count));
                fct_chk_eq_int(count, ncomponents);
            }

            // Close the file
            fct_req(0 == esio_file_close(state));

            // Free the temporary
            free(value);

            // Reopen the file using normal HDF5 APIs on root processor
            // Examine the contents to ensure it matches
            if (world_rank == 0) {
                value = calloc(ncomponents, sizeof(TYPE));
                fct_req(value);
                const hid_t file_id
                    = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
                fct_req(0 <= AFFIX(H5LTget_attribute)(file_id, "/",
                                                      "attribute", value));
                for (int i = 0; i < ncomponents; ++i) {
                    fct_chk(i + 5678 == value[i]);
                }
                fct_req(0 <= H5Fclose(file_id));
                free(value);
            }

            // Re-read the file in a distributed manner and verify contents
            value = calloc(ncomponents, sizeof(TYPE));
            fct_req(value);
            fct_req(0 == esio_file_open(state, filename, 0));
            fct_req(0 == AFFIX(esio_attribute_readv)(state, "attribute",
                                                     value, ncomponents));
            for (int i = 0; i < ncomponents; ++i) {
                fct_chk(i + 5678 == value[i]);
            }
            free(value);
            fct_req(0 == esio_file_close(state));
        }
        FCT_TEST_END();

    }
    FCT_FIXTURE_SUITE_END();
}
FCT_END()
