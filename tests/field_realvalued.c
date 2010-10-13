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

#include "fct.h"

FCT_BGN()
{
    // MPI setup: MPI_Init and atexit(MPI_Finalize)
    int world_size, world_rank;
    MPI_Init(&argc, &argv);
    atexit((void (*) ()) MPI_Finalize);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

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

            // Only one process generates a unique filename
            int filenamelen;
            if (world_rank == 0) {
                filename = tempnam(output_dir, "l1tst");
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
            int na, ast, asz;
            int nb, bst, bsz;
            int nc, cst, csz;

            for (int casenum = 0; casenum < 3; ++casenum) {

                ESIO_MPICHKQ(MPI_Barrier(MPI_COMM_WORLD)); // Barrier

                na = ast = asz = nb = bst = bsz = nc = cst = csz = -1;
                switch (casenum) {
                    case 0:
                        // Data uniformly partitioned in the fastest index
                        asz = 3;
                        na  = asz * world_size;
                        ast = asz * world_rank;

                        nb  = bsz = 2;
                        nc  = csz = 2;
                        bst = cst = 0;
                        break;
                    case 1:
                        // Data uniformly partitioned in the medium index
                        bsz = 3;
                        nb  = bsz * world_size;
                        bst = bsz * world_rank;

                        na  = asz = 2;
                        nc  = csz = 2;
                        ast = cst = 0;
                        break;
                    case 2:
                        // Data uniformly partitioned in the slow index
                        csz = 3;
                        nc  = csz * world_size;
                        cst = csz * world_rank;

                        na  = asz = 2;
                        nb  = bsz = 2;
                        ast = bst = 0;
                        break;
                    default:
                        fct_req(0); // Sanity failure
                }

                fct_req(0 == esio_file_create(state, filename, 1));

                const size_t nelem = asz * bsz * csz;
                TEST_REAL * field  = calloc(nelem, sizeof(TEST_REAL));
                fct_req(field);

                { // Populate the field with test data
                    TEST_REAL * p_field = field;
                    for (int k = cst; k < cst + csz; ++k) {
                        for (int j = bst; j < bst + bsz; ++j) {
                            for (int i = ast; i < ast + asz; ++i) {
                                *p_field++ = 2*(i+3)+5*(j+7)+11*(k+13);
                            }
                        }
                    }
                }

                // Write the data to disk
                const int status = TEST_ESIO_FIELD_WRITE(state, "field", field,
                        na, ast, asz, nb, bst, bsz, nc, cst, csz);
                fct_req(status == 0);

                // Close the file
                fct_req(0 == esio_file_close(state));

                // Free the field
                free(field);

                // Reopen the file using normal HDF5 APIs on root processor
                // Examine the contents to ensure it matches
                if (world_rank == 0) {
                    field = calloc(na * nb * nc, sizeof(TEST_REAL));
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
                                    value, 2*(i+3)+5*(j+7)+11*(k+13));
                            }
                        }
                    }

                    const herr_t status2 = H5Fclose(file_id);
                    fct_req(status2 >= 0);
                    free(field);

                    const int status = unlink(filename);
                    assert(0 == status);
                }
            }

        }
        FCT_TEST_END();

    }
    FCT_FIXTURE_SUITE_END();
}
FCT_END();
