//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// esio 0.1.1: ExaScale IO library for turbulence simulation restart files
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

#if    !defined(REAL) \
    || !defined(REAL_H5T) \
    || !defined(AFFIX)
#error "line_template.c should be #included after appropriate #defines"
#endif

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
        "--partitioned-size",
        "-p",
        FCTCL_STORE_VALUE,
        "Sets size/rank of directions to be uniformly partitioned"
    },
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
    {
        "--auxstride-a",
        NULL,
        FCTCL_STORE_VALUE,
        "Auxiliary stride added to the a direction"
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
    const int partitioned_size = (int) strtol(
        fctcl_val2("--partitioned-size","3"), (char **) NULL, 10);
    if (partitioned_size < 1) {
        fprintf(stderr, "\n--partitioned-size=%d < 1\n", partitioned_size);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    const int ncomponents = (int) strtol(
        fctcl_val2("--ncomponents","2"), (char **) NULL, 10);
    if (ncomponents < 1) {
        fprintf(stderr, "\n--ncomponents=%d < 1\n", ncomponents);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Data uniformly partitioned in the fastest (only) index
    const int alocal  = partitioned_size;
    const int aglobal = alocal * world_size;
    const int astart  = alocal * world_rank;

    // Retrieve and process auxilary stride options
    const int auxstride_a = (int) strtol(
        fctcl_val2("--auxstride-a","0"), (char **) NULL, 10);
    if (auxstride_a < 0) {
        fprintf(stderr, "\n--auxstride-a=%d < 0\n", auxstride_a);
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
    char * filetemplate = create_testfiletemplate(output_dir, __FILE__);
    (void) input_dir;  // Possibly unused
    char * filename = NULL;
    esio_handle handle;

    FCT_FIXTURE_SUITE_BGN(line_suite)
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

            // Initialize ESIO handle
            handle = esio_handle_initialize(MPI_COMM_WORLD);
            assert(handle);
        }
        FCT_SETUP_END();

        FCT_TEARDOWN_BGN()
        {
            // Finalize ESIO handle
            esio_handle_finalize(handle);

            // Clean up the unique file and filename
            if (world_rank == 0) {
                if (!preserve) unlink(filename);
            }
            if (filename) free(filename);
        }
        FCT_TEARDOWN_END();

        // Test scalar-valued lines, including overwrite details
        FCT_TEST_BGN(line)
        {
            REAL *line;

            // Compute stride information for each direction;
            // strides expressed in sizeof(REAL) which happens to be 1.
            const int astride  = 1 + auxstride_a * 1;
            const size_t nelem = alocal * astride;

            // Allocate storage for local portion of global line
            line = calloc(nelem, sizeof(REAL));
            fct_req(line);

            // Establish and then check the parallel decomposition
            fct_req(0 == esio_line_establish(handle, aglobal, astart, alocal));
            {
                int tmp_aglobal, tmp_astart, tmp_alocal;
                fct_req(0 == esio_line_established(handle, &tmp_aglobal,
                                                           &tmp_astart,
                                                           &tmp_alocal));
                fct_chk_eq_int(aglobal, tmp_aglobal);
                fct_chk_eq_int(astart,  tmp_astart);
                fct_chk_eq_int(alocal,  tmp_alocal);

                fct_req(0 == esio_line_established(handle, NULL, NULL, NULL));
            }

            // Open file
            fct_req(0 == esio_file_create(handle, filename, 1));

            // Write zeros to disk and flush the buffers
            fct_req(0 == AFFIX(esio_line_write)(
                               handle, "line", line, astride));
            fct_req(0 == esio_file_flush(handle));

            // Populate local line with test data
            for (int i = 0; i < alocal; ++i) {
                const REAL value = (REAL) 2*((i + astart) +  3);
                line[i*astride] = value;
            }

            // Overwrite zeros on disk with test data
            if (auxstride_a) {
                // Noncontiguous; exercise non-default stride arguments
                fct_req(0 == AFFIX(esio_line_write)(
                                   handle, "line", line, astride));
            } else {
                // Contiguous; exercise default stride arguments
                fct_req(0 == AFFIX(esio_line_write)(
                                   handle, "line", line, 0));
            }

            { // Ensure the global size was written correctly
                int tmp_aglobal;
                fct_req(0 == esio_line_size(handle, "line", &tmp_aglobal));
                fct_chk_eq_int(aglobal, tmp_aglobal);
            }

            // Close the file
            fct_req(0 == esio_file_close(handle));

            // Free the line
            free(line);

            // Reopen the file using normal HDF5 APIs on root processor
            // Examine the contents (contiguously) to ensure it matches
            if (world_rank == 0) {
                line = calloc(aglobal, sizeof(REAL));
                fct_req(line);
                const hid_t file_id
                    = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
                fct_req(0 <= H5LTread_dataset(file_id, "line",
                                              REAL_H5T, line));

                REAL *p_line = line;
                for (int i = 0; i < aglobal; ++i) {
                    const REAL expected = (REAL) 2*(i+3);
                    const REAL value    = *p_line++;
                    fct_chk_eq_dbl(value, expected);
                }

                fct_req(0 <= H5Fclose(file_id));
                free(line);
            }

            // Re-read the file in a distributed manner and verify contents
            // TODO Ensure non-referenced memory locations remain unmodified
            line = calloc(nelem, sizeof(REAL));
            fct_req(line);
            fct_req(0 == esio_file_open(handle, filename, 0));
            if (auxstride_a) {
                // Noncontiguous; exercise non-default stride arguments
                fct_req(0 == AFFIX(esio_line_read)(
                                   handle, "line", line, astride));
            } else {
                // Contiguous; exercise default stride arguments
                fct_req(0 == AFFIX(esio_line_read)(
                                   handle, "line", line, 0));
            }
            for (int i = 0; i < alocal; ++i) {
                const REAL expected = (REAL) 2*((i + astart) +  3);
                const REAL value = line[i*astride];
                fct_chk_eq_dbl(value, expected);
            }
            free(line);
            fct_req(0 == esio_file_close(handle));
        }
        FCT_TEST_END();

        // Test vector-valued lines
        FCT_TEST_BGN(vline)
        {
            REAL *vline;

            // Compute stride information for each direction;
            // strides expressed in sizeof(REAL).
            const int astride  = ncomponents + auxstride_a * ncomponents;
            const size_t nelem = alocal * astride;

            // Allocate storage for local portion of global vline
            vline = calloc(nelem, sizeof(REAL));
            fct_req(vline);

            // Establish the parallel decomposition
            fct_req(0 == esio_line_establish(handle, aglobal, astart, alocal));

            // Open file
            fct_req(0 == esio_file_create(handle, filename, 1));

            // Populate local vline with test data
            for (int i = 0; i < alocal; ++i) {
                for (int h = 0; h < ncomponents; ++h) {
                    vline[i*astride + h] = (REAL) 2*((i + astart) +  3) - h;
                }
            }

            // Write vline to file
            if (auxstride_a) {
                // Noncontiguous; exercise non-default stride arguments
                fct_req(0 == AFFIX(esio_line_writev)(
                                   handle, "vline", vline,
                                   astride, ncomponents));
            } else {
                // Contiguous; exercise default stride arguments
                fct_req(0 == AFFIX(esio_line_writev)(
                                   handle, "vline", vline,
                                   0, ncomponents));
            }

            { // Ensure the global size was written correctly
                int tmp_aglobal, tmp_ncomponents;
                fct_req(0 == esio_line_sizev(handle, "vline",
                                             &tmp_aglobal,
                                             &tmp_ncomponents));

                fct_chk_eq_int(aglobal,     tmp_aglobal);
                fct_chk_eq_int(ncomponents, tmp_ncomponents);
            }

            // Close the file
            fct_req(0 == esio_file_close(handle));

            // Free the vline
            free(vline);

            // Reopen the file using normal HDF5 APIs on root processor
            // Examine the contents to ensure it matches
            if (world_rank == 0) {
                vline = calloc(aglobal * ncomponents, sizeof(REAL));
                fct_req(vline);
                const hid_t file_id
                    = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
                hid_t type_id = REAL_H5T;
                if (ncomponents > 1) {
                    const hsize_t dims[1] = { ncomponents };
                    type_id = H5Tarray_create2(REAL_H5T, 1, dims);
                    fct_req(type_id >= 0);
                }
                fct_req(0 <= H5LTread_dataset(file_id, "vline",
                                              type_id, vline));
                if (ncomponents > 1) H5Tclose(type_id);

                REAL *p_line = vline;
                for (int i = 0; i < aglobal; ++i) {
                    for (int h = 0; h < ncomponents; ++h) {
                        const REAL value = *p_line++;
                        fct_chk_eq_dbl(value, (REAL) 2*(i+3)-h);
                    }
                }

                fct_req(0 <= H5Fclose(file_id));
                free(vline);
            }

            // Re-read the file in a distributed manner and verify contents
            // TODO Ensure non-referenced memory locations remain unmodified
            vline = calloc(nelem, sizeof(REAL));
            fct_req(vline);
            fct_req(0 == esio_file_open(handle, filename, 0));
            if (auxstride_a) {
                // Noncontiguous; exercise non-default stride arguments
                fct_req(0 == AFFIX(esio_line_readv)(
                                   handle, "vline", vline,
                                   astride, ncomponents));
            } else {
                // Contiguous; exercise default stride arguments
                fct_req(0 == AFFIX(esio_line_readv)(
                                   handle, "vline", vline,
                                   0, ncomponents));
            }
            for (int i = 0; i < alocal; ++i) {
                for (int h = 0; h < ncomponents; ++h) {
                    const REAL expected = (REAL) 2*((i + astart) +  3) - h;
                    const REAL value = vline[i*astride + h];
                    fct_chk_eq_dbl(value, expected);
                }
            }
            free(vline);
            fct_req(0 == esio_file_close(handle));
        }
        FCT_TEST_END();

    }
    FCT_FIXTURE_SUITE_END();

    if (filetemplate) free(filetemplate);
}
FCT_END()
