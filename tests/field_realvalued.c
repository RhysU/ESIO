/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PSDNS Development Team
 *
 * Please see https://wiki.ices.utexas.edu/PSDNS for more information.
 *
 * This file is part of the esio library.
 *
 * ESIO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESIO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ESIO.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * field_realvalued.c: real-valued unit test template for ESIO
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef TEST_REAL
// You probably wanted to compile field_double.c or field_float.c
#error Need to #define TEST_REAL prior to compiling field_realvalued.c
#endif

// See 'feature_test_macros(7) for details'
#define _GNU_SOURCE

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
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
static fctcl_init_t my_cl_options[] = {
    {
        "--partitioned-size",
        "-p",
        FCTCL_STORE_VALUE,
        "Sets size/rank of directions to be uniformly partitioned"
    },
    {
        "--unpartitioned-size",
        "-u",
        FCTCL_STORE_VALUE,
        "Sets size of directions which are not partitioned"
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

    // Retrieve and sanity check command line options
    preserve = fctcl_is("--preserve");
    const int partitioned_size = (int) strtol(
        fctcl_val2("--partitioned-size","3"), (char **) NULL, 10);
    if (partitioned_size < 1) {
        fprintf(stderr, "\n--partitioned-size=%d < 1\n", partitioned_size);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    const int unpartitioned_size = (int) strtol(
        fctcl_val2("--unpartitioned-size","2"), (char **) NULL, 10);
    if (unpartitioned_size < 1) {
        fprintf(stderr, "\n--unpartitioned-size=%d < 1\n", unpartitioned_size);
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
    char * filename = NULL;
    esio_state state;

    FCT_FIXTURE_SUITE_BGN(field_realvalued)
    {
        FCT_SETUP_BGN()
        {
            (void) input_dir; // Unused
            H5Eset_auto2(H5E_DEFAULT, hdf5_handler, hdf5_client_data);
            esio_set_error_handler(esio_handler);

            state = esio_init(MPI_COMM_WORLD);
            assert(state);
        }
        FCT_SETUP_END();

        FCT_TEARDOWN_BGN()
        {
            esio_finalize(state);
            ESIO_MPICHKQ(MPI_Barrier(MPI_COMM_WORLD)); // Barrier
        }
        FCT_TEARDOWN_END();

        FCT_TEST_BGN(uniform_directionally_split)
        {
            int nc, cst, csz;
            int nb, bst, bsz;
            int na, ast, asz;

            for (int casenum = 0; casenum < 3; ++casenum) {

                ESIO_MPICHKQ(MPI_Barrier(MPI_COMM_WORLD)); // Barrier

                na = ast = asz = nb = bst = bsz = nc = cst = csz = -1;
                switch (casenum) {
                    case 0: // Data uniformly partitioned in the fastest index
                        asz = partitioned_size;
                        na  = asz * world_size;
                        ast = asz * world_rank;

                        nb  = bsz = unpartitioned_size;
                        nc  = csz = unpartitioned_size;
                        bst = cst = 0;
                        break;

                    case 1: // Data uniformly partitioned in the medium index
                        bsz = partitioned_size;
                        nb  = bsz * world_size;
                        bst = bsz * world_rank;

                        na  = asz = unpartitioned_size;
                        nc  = csz = unpartitioned_size;
                        ast = cst = 0;
                        break;

                    case 2: // Data uniformly partitioned in the slow index
                        csz = partitioned_size;
                        nc  = csz * world_size;
                        cst = csz * world_rank;

                        na  = asz = unpartitioned_size;
                        nb  = bsz = unpartitioned_size;
                        ast = bst = 0;
                        break;
                    default:
                        fct_req(0); // Sanity failure
                }

                // Rank 0 generates a unique filename and broadcasts it
                int filenamelen;
                if (world_rank == 0) {
                    filename = tempnam(output_dir, "l00t");
                    if (preserve) {
                        printf("\ncase %d filename: %s\n", casenum, filename);
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

                // Open unique file
                fct_req(0 == esio_file_create(state, filename, 1));

                // Determine local rank's portion of global field
                const size_t nelem = asz * bsz * csz;
                TEST_REAL * field  = calloc(nelem, sizeof(TEST_REAL));
                fct_req(field);

                { // Populate local field with test data
                    TEST_REAL * p_field = field;
                    for (int k = cst; k < cst + csz; ++k) {
                        for (int j = bst; j < bst + bsz; ++j) {
                            for (int i = ast; i < ast + asz; ++i) {
                                *p_field++
                                    = (TEST_REAL) 2*(i+3)+5*(j+7)+11*(k+13);
                            }
                        }
                    }
                }

                // Write local data to disk
                int status = TEST_ESIO_FIELD_WRITE(state, "field", field,
                        nc, cst, csz, nb, bst, bsz, na, ast, asz);
                fct_req(status == 0);

                // Ensure the global size was written correctly
                {
                    int tmp_na, tmp_nb, tmp_nc;
                    fct_req(0 == esio_field_size(state, "field",
                                                 &tmp_nc, &tmp_nb, &tmp_na));
                    fct_chk_eq_int(nc, tmp_nc);
                    fct_chk_eq_int(nb, tmp_nb);
                    fct_chk_eq_int(na, tmp_na);
                }

                // Close the file
                fct_req(0 == esio_file_close(state));

                // Free the field
                free(field);

                // Reopen the file using normal HDF5 APIs on root processor
                // Examine the contents to ensure it matches
                if (world_rank == 0) {
                    field = calloc(nc * nb * na, sizeof(TEST_REAL));
                    fct_req(field);
                    const hid_t file_id
                        = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
                    const herr_t status1
                        = TEST_H5LTREAD_DATASET(file_id, "field", field);
                    fct_req(status1 >= 0);

                    TEST_REAL * p_field = field;
                    for (int k = 0; k < nc; ++k) {
                        for (int j = 0; j < nb; ++j) {
                            for (int i = 0; i < na; ++i) {
                                const TEST_REAL value = *p_field++;
                                fct_chk_eq_dbl(
                                        value,
                                        (TEST_REAL) 2*(i+3)+5*(j+7)+11*(k+13));
                            }
                        }
                    }

                    const herr_t status2 = H5Fclose(file_id);
                    fct_req(status2 >= 0);
                    free(field);
                }

                // Re-read the file in a distributed manner and verify contents
                field = calloc(nelem, sizeof(TEST_REAL));
                fct_req(field);
                fct_req(0 == esio_file_open(state, filename, 0));
                status = TEST_ESIO_FIELD_READ(state, "field", field,
                        nc, cst, csz, nb, bst, bsz, na, ast, asz);
                fct_req(status == 0);
                {
                    TEST_REAL * p_field = field;
                    for (int k = cst; k < cst + csz; ++k) {
                        for (int j = bst; j < bst + bsz; ++j) {
                            for (int i = ast; i < ast + asz; ++i) {
                                const TEST_REAL value = *p_field++;
                                fct_chk_eq_dbl(
                                    value,
                                    (TEST_REAL) 2*(i+3)+5*(j+7)+11*(k+13));
                            }
                        }
                    }
                }
                free(field);
                fct_req(0 == esio_file_close(state));

                // Clean up
                if (world_rank == 0) {
                    if (!preserve) fct_req(0 == unlink(filename));
                }
                if (filename) free(filename);
            }

        }
        FCT_TEST_END();

    }
    FCT_FIXTURE_SUITE_END();
}
FCT_END()
