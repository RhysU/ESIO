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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <mpi.h>
#include <hdf5.h>
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

    // Retrieve options
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
    char * filename = NULL;
    esio_handle state;

    FCT_FIXTURE_SUITE_BGN(esio_file)
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

        FCT_TEST_BGN(success_code)
        {
            fct_chk_eq_int(0, ESIO_SUCCESS); // Success is zero
            fct_chk(!ESIO_SUCCESS);          // Not success is true
        }
        FCT_TEST_END();

        FCT_TEST_BGN(file_create_and_open)
        {
            // No open file so esio_file_path is empty
            fct_req(NULL == esio_file_path(state));

            // Create with overwrite should always work
            fct_req(0 == esio_file_create(state, filename, 1 /* overwrite */));

            // Flush flush flush should always work
            fct_req(0 == esio_file_flush(state));
            fct_req(0 == esio_file_flush(state));
            fct_req(0 == esio_file_flush(state));

            // Check that esio_file_path points to a valid canonical location
            struct stat statbuf;
            char *file_path = esio_file_path(state);
            fct_req(file_path != NULL);
            fct_chk_startswith_str(file_path, "/"); // Absolute?
            fct_req(0 == stat(file_path, &statbuf));
            free(file_path);

            // Close the file
            fct_req(0 == esio_file_close(state));

            // Double closure should silently succeed
            fct_req(0 == esio_file_close(state));

            // No open file so esio_file_path is empty
            fct_req(NULL == esio_file_path(state));

            // Create without overwrite should fail
            H5Eset_auto(H5E_DEFAULT, NULL, NULL);
            esio_set_error_handler_off();
            fct_req(0 != esio_file_create(state, filename, 0 /* no overwrite */));
            esio_set_error_handler(esio_handler);
            H5Eset_auto2(H5E_DEFAULT, hdf5_handler, hdf5_client_data);

            // Remove the file and create without overwrite should succeed
            ESIO_MPICHKQ(MPI_Barrier(MPI_COMM_WORLD)); // Synchronize
            if (world_rank == 0) {
                fct_req(0 == unlink(filename));
            }
            ESIO_MPICHKQ(MPI_Barrier(MPI_COMM_WORLD)); // Synchronize
            fct_req(0 == esio_file_create(state, filename, 0 /* no overwrite */));

            // Close the file
            fct_req(0 == esio_file_close(state));

            // Ensure we can open in read-only mode
            fct_req(0 == esio_file_open(state, filename, 0 /* read-only */));

            // Close the file
            fct_req(0 == esio_file_close(state));

            // Ensure we can open in read-write mode
            fct_req(0 == esio_file_open(state, filename, 1 /* read-write */));

            // Close the file
            fct_req(0 == esio_file_close(state));
        }
        FCT_TEST_END();

        FCT_TEST_BGN(file_clone)
        {
            // Dynamically create the empty.h5 source filename per input_dir
            if (input_dir == NULL && world_rank == 0) {
                fprintf(stderr, "\nESIO_TEST_INPUT_DIR not in environment!\n");
            }
            fct_req(NULL != input_dir /* check ESIO_TEST_INPUT_DIR set */);
            int srcfilelen = strlen(input_dir);
            fct_req(srcfilelen > 0);
            srcfilelen += strlen("/empty.h5");
            srcfilelen += 1;
            char * srcfile = calloc(srcfilelen, 1);
            fct_req(NULL != srcfile);
            fct_req(NULL != strcat(srcfile, input_dir));
            fct_req(NULL != strcat(srcfile, "/empty.h5"));

            // Clone with overwrite should always work
            fct_req(0 == esio_file_clone(state, srcfile, filename, 1));

            // Close the file
            fct_req(0 == esio_file_close(state));

            // Create without overwrite should fail
            H5Eset_auto(H5E_DEFAULT, NULL, NULL);
            esio_set_error_handler_off();
            fct_req(0 != esio_file_clone(state, srcfile, filename, 0));
            esio_set_error_handler(esio_handler);
            H5Eset_auto2(H5E_DEFAULT, hdf5_handler, hdf5_client_data);

            // Remove the file and create without overwrite should succeed
            ESIO_MPICHKQ(MPI_Barrier(MPI_COMM_WORLD)); // Synchronize
            if (world_rank == 0) {
                fct_req(0 == unlink(filename));
            }
            ESIO_MPICHKQ(MPI_Barrier(MPI_COMM_WORLD)); // Synchronize
            fct_req(0 == esio_file_clone(state, srcfile, filename, 0));

            // Close the file
            fct_req(0 == esio_file_close(state));

            // Deallocate the srcfile name
            if (srcfile) free(srcfile);
        }
        FCT_TEST_END();

        // Much of this functionality is covered in restart_helpers.c
        // and restart_rename.{sh,c}.  This covers only the public API's
        // most basic usage.
        FCT_TEST_BGN(file_close_restart)
        {
            struct stat statbuf;

            // Form template and expected file paths from temporary filename
            char *template = malloc(strlen(filename) + 2);
            fct_req(template);
            strcpy(template, filename);
            strcat(template, "#");
            char *restart0 = malloc(strlen(filename) + 2);
            fct_req(restart0);
            strcpy(restart0, filename);
            strcat(restart0, "0");
            char *restart1 = malloc(strlen(filename) + 2);
            fct_req(restart1);
            strcpy(restart1, filename);
            strcat(restart1, "1");

            ESIO_MPICHKQ(MPI_Barrier(MPI_COMM_WORLD)); // Synchronize

            // No files should exist prior to the test kicking off
            fct_req(-1 == stat(filename, &statbuf));
            fct_req(-1 == stat(restart0, &statbuf));
            fct_req(-1 == stat(restart1, &statbuf));

            ESIO_MPICHKQ(MPI_Barrier(MPI_COMM_WORLD)); // Synchronize

            // Create with overwrite should always work
            fct_req( 0 == esio_file_create(state, filename, 1 /* overwrite */));
            fct_req( 0 == esio_file_close_restart(state, template, 2));
            fct_req(-1 == stat(filename, &statbuf));
            fct_req( 0 == stat(restart0, &statbuf));
            fct_req(-1 == stat(restart1, &statbuf));

            ESIO_MPICHKQ(MPI_Barrier(MPI_COMM_WORLD)); // Synchronize

            // Create without overwrite should work now
            fct_req(0 == esio_file_create(state, filename, 0));
            fct_req(0 == esio_file_close_restart(state, template, 2));
            fct_req(-1 == stat(filename, &statbuf));
            fct_req( 0 == stat(restart0, &statbuf));
            fct_req( 0 == stat(restart1, &statbuf));

            ESIO_MPICHKQ(MPI_Barrier(MPI_COMM_WORLD)); // Synchronize

            // Clean up
            if (world_rank == 0) {
                fct_req(0 == unlink(restart0));
                fct_req(0 == unlink(restart1));
            }
            free(template);
            free(restart0);
            free(restart1);
        }
        FCT_TEST_END();

    }
    FCT_FIXTURE_SUITE_END();

    if (filetemplate) free(filetemplate);
}
FCT_END()
