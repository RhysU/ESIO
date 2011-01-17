//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// esio 0.1.2: ExaScale IO library for turbulence simulation restart files
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
#error "plane_template.c should be #included after appropriate #defines"
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

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
        "--distribute",
        "-d",
        FCTCL_STORE_VALUE,
        "Choose which direction (b or a) is distributed"
    },
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
        "--ncomponents",
        "-n",
        FCTCL_STORE_VALUE,
        "Sets the number of components used in vplane test"
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
    {
        "--auxstride-b",
        NULL,
        FCTCL_STORE_VALUE,
        "Auxiliary stride added to the b direction"
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
    const int unpartitioned_size = (int) strtol(
        fctcl_val2("--unpartitioned-size","2"), (char **) NULL, 10);
    if (unpartitioned_size < 1) {
        fprintf(stderr, "\n--unpartitioned-size=%d < 1\n", unpartitioned_size);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    const int ncomponents = (int) strtol(
        fctcl_val2("--ncomponents","2"), (char **) NULL, 10);
    if (ncomponents < 1) {
        fprintf(stderr, "\n--ncomponents=%d < 1\n", ncomponents);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Retrieve and process parallel distribution options
    int bglobal, bstart, blocal;
    int aglobal, astart, alocal;
    const char *distribute = fctcl_val2("--distribute", "a");
    if (fctstr_ieq(distribute,"b")) {
        // Data uniformly partitioned in the medium index
        blocal = partitioned_size;
        bglobal = blocal * world_size;
        bstart = blocal * world_rank;

        aglobal = alocal = unpartitioned_size;
        astart = 0;
    } else if (fctstr_ieq(distribute,"a")) {
        // Data uniformly partitioned in the fastest index
        alocal = partitioned_size;
        aglobal = alocal * world_size;
        astart = alocal * world_rank;

        bglobal = blocal = unpartitioned_size;
        bstart = 0;
    } else {
        fprintf(stderr, "\n--distribute=%s not recognized\n", distribute);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Retrieve and process auxilary stride options
    const int auxstride_a = (int) strtol(
        fctcl_val2("--auxstride-a","0"), (char **) NULL, 10);
    if (auxstride_a < 0) {
        fprintf(stderr, "\n--auxstride-a=%d < 0\n", auxstride_a);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    const int auxstride_b = (int) strtol(
        fctcl_val2("--auxstride-b","0"), (char **) NULL, 10);
    if (auxstride_b < 0) {
        fprintf(stderr, "\n--auxstride-b=%d < 0\n", auxstride_b);
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

    FCT_FIXTURE_SUITE_BGN(plane_suite)
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

            // Initialize ESIO handle
            handle = esio_handle_initialize(MPI_COMM_WORLD);
            fct_req(handle);
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

        // Test scalar-valued planes, including overwrite details
        FCT_TEST_BGN(plane)
        {
            REAL *plane;

            // Compute stride information for each direction;
            // strides expressed in sizeof(REAL) which happens to be 1.
            const int astride  = 1                + auxstride_a * 1;
            const int bstride  = alocal * astride + auxstride_b * 1;
            const size_t nelem = blocal * bstride;

            // Allocate storage for local portion of global plane
            plane = calloc(nelem, sizeof(REAL));
            fct_req(plane);

            // Establish and then check the parallel decomposition
            fct_req(0 == esio_plane_establish(handle, bglobal, bstart, blocal,
                                                      aglobal, astart, alocal));
            {
                int tmp_bglobal, tmp_bstart, tmp_blocal;
                int tmp_aglobal, tmp_astart, tmp_alocal;
                fct_req(0 == esio_plane_established(handle, &tmp_bglobal,
                                                            &tmp_bstart,
                                                            &tmp_blocal,
                                                            &tmp_aglobal,
                                                            &tmp_astart,
                                                            &tmp_alocal));
                fct_chk_eq_int(bglobal, tmp_bglobal);
                fct_chk_eq_int(bstart,  tmp_bstart);
                fct_chk_eq_int(blocal,  tmp_blocal);
                fct_chk_eq_int(aglobal, tmp_aglobal);
                fct_chk_eq_int(astart,  tmp_astart);
                fct_chk_eq_int(alocal,  tmp_alocal);

                fct_req(0 == esio_plane_established(handle, NULL, NULL, NULL,
                                                            NULL, NULL, NULL));
            }

            // Open file
            fct_req(0 == esio_file_create(handle, filename, 1));

            // Write zeros to disk and flush the buffers
            fct_req(0 == AFFIX(esio_plane_write)(
                                handle, "plane", plane, bstride, astride));
            fct_req(0 == esio_file_flush(handle));

            // Populate local plane with test data
            for (int j = 0; j < blocal; ++j) {
                for (int i = 0; i < alocal; ++i) {
                    const REAL value =  (REAL) 5*((j + bstart) +  7)
                                            +  2*((i + astart) +  3);
                    plane[j*bstride + i*astride] = value;
                }
            }

            // Overwrite zeros on disk with test data
            if (auxstride_a || auxstride_b) {
                // Noncontiguous; exercise non-default stride arguments
                fct_req(0 == AFFIX(esio_plane_write)(
                                    handle, "plane", plane, bstride, astride));
            } else {
                // Contiguous; exercise default stride arguments
                fct_req(0 == AFFIX(esio_plane_write)(
                                    handle, "plane", plane, 0, 0));
            }

            { // Ensure the global size was written correctly
                int tmp_aglobal, tmp_bglobal;
                fct_req(0 == esio_plane_size(handle, "plane",
                                             &tmp_bglobal,
                                             &tmp_aglobal));
                fct_chk_eq_int(bglobal, tmp_bglobal);
                fct_chk_eq_int(aglobal, tmp_aglobal);
            }

            // Close the file
            fct_req(0 == esio_file_close(handle));

            // Free the plane
            free(plane);

            // Reopen the file using normal HDF5 APIs on root processor
            // Examine the contents (contiguously) to ensure it matches
            if (world_rank == 0) {
                plane = calloc(bglobal * aglobal, sizeof(REAL));
                fct_req(plane);
                const hid_t file_id
                    = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
                fct_req(0 <= H5LTread_dataset(file_id, "plane",
                                              REAL_H5T, plane));

                REAL *p_plane = plane;
                for (int j = 0; j < bglobal; ++j) {
                    for (int i = 0; i < aglobal; ++i) {
                        const REAL expected = (REAL) 2*(i+3)+5*(j+7);
                        const REAL value    = *p_plane++;
                        fct_chk_eq_dbl(value, expected);
                    }
                }

                fct_req(0 <= H5Fclose(file_id));
                free(plane);
            }

            // Re-read the file in a distributed manner and verify contents
            // TODO Ensure non-referenced memory locations remain unmodified
            plane = calloc(nelem, sizeof(REAL));
            fct_req(plane);
            fct_req(0 == esio_file_open(handle, filename, 0));
            if (auxstride_a || auxstride_b) {
                // Noncontiguous; exercise non-default stride arguments
                fct_req(0 == AFFIX(esio_plane_read)(
                                    handle, "plane", plane, bstride, astride));
            } else {
                // Contiguous; exercise default stride arguments
                fct_req(0 == AFFIX(esio_plane_read)(
                                    handle, "plane", plane, 0, 0));
            }
            for (int j = 0; j < blocal; ++j) {
                for (int i = 0; i < alocal; ++i) {
                    const REAL expected =  (REAL) 5*((j + bstart) +  7)
                                               +  2*((i + astart) +  3);
                    const REAL value = plane[j*bstride + i*astride];
                    fct_chk_eq_dbl(value, expected);
                }
            }
            free(plane);
            fct_req(0 == esio_file_close(handle));
        }
        FCT_TEST_END();

        // Test vector-valued planes
        FCT_TEST_BGN(vplane)
        {
            REAL *vplane;

            // Compute stride information for each direction;
            // strides expressed in sizeof(REAL).
            const int astride  = ncomponents      + auxstride_a * ncomponents;
            const int bstride  = alocal * astride + auxstride_b * ncomponents;
            const size_t nelem = blocal * bstride;

            // Allocate storage for local portion of global vplane
            vplane = calloc(nelem, sizeof(REAL));
            fct_req(vplane);

            // Establish the parallel decomposition
            fct_req(0 == esio_plane_establish(handle, bglobal, bstart, blocal,
                                                      aglobal, astart, alocal));

            // Open file
            fct_req(0 == esio_file_create(handle, filename, 1));

            // Populate local vplane with test data
            for (int j = 0; j < blocal; ++j) {
                for (int i = 0; i < alocal; ++i) {
                    for (int h = 0; h < ncomponents; ++h) {
                        vplane[j*bstride + i*astride + h]
                                =  (REAL) 5*((j + bstart) +  7)
                                       +  2*((i + astart) +  3)
                                       - h;
                    }
                }
            }

            // Write vplane to file
            if (auxstride_a || auxstride_b) {
                // Noncontiguous; exercise non-default stride arguments
                fct_req(0 == AFFIX(esio_plane_writev)(
                                    handle, "vplane", vplane,
                                    bstride, astride, ncomponents));
            } else {
                // Contiguous; exercise default stride arguments
                fct_req(0 == AFFIX(esio_plane_writev)(
                                    handle, "vplane", vplane,
                                    0, 0, ncomponents));
            }

            { // Ensure the global size was written correctly
                int tmp_aglobal, tmp_bglobal, tmp_ncomponents;
                fct_req(0 == esio_plane_sizev(handle, "vplane",
                                              &tmp_bglobal,
                                              &tmp_aglobal,
                                              &tmp_ncomponents));

                fct_chk_eq_int(bglobal,     tmp_bglobal);
                fct_chk_eq_int(aglobal,     tmp_aglobal);
                fct_chk_eq_int(ncomponents, tmp_ncomponents);
            }

            // Close the file
            fct_req(0 == esio_file_close(handle));

            // Free the vplane
            free(vplane);

            // Reopen the file using normal HDF5 APIs on root processor
            // Examine the contents to ensure it matches
            if (world_rank == 0) {
                vplane = calloc(bglobal * aglobal * ncomponents,
                                sizeof(REAL));
                fct_req(vplane);
                const hid_t file_id
                    = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
                hid_t type_id = REAL_H5T;
                if (ncomponents > 1) {
                    const hsize_t dims[1] = { ncomponents };
                    type_id = H5Tarray_create2(REAL_H5T, 1, dims);
                    fct_req(type_id >= 0);
                }
                fct_req(0 <= H5LTread_dataset(file_id, "vplane",
                                              type_id, vplane));
                if (ncomponents > 1) H5Tclose(type_id);

                REAL *p_plane = vplane;
                for (int j = 0; j < bglobal; ++j) {
                    for (int i = 0; i < aglobal; ++i) {
                        for (int h = 0; h < ncomponents; ++h) {
                            const REAL value = *p_plane++;
                            fct_chk_eq_dbl(
                                value,
                                (REAL) 2*(i+3)+5*(j+7)-h);
                        }
                    }
                }

                fct_req(0 <= H5Fclose(file_id));
                free(vplane);
            }

            // Re-read the file in a distributed manner and verify contents
            // TODO Ensure non-referenced memory locations remain unmodified
            vplane = calloc(nelem, sizeof(REAL));
            fct_req(vplane);
            fct_req(0 == esio_file_open(handle, filename, 0));
            if (auxstride_a || auxstride_b) {
                // Noncontiguous; exercise non-default stride arguments
                fct_req(0 == AFFIX(esio_plane_readv)(
                                    handle, "vplane", vplane,
                                    bstride, astride, ncomponents));
            } else {
                // Contiguous; exercise default stride arguments
                fct_req(0 == AFFIX(esio_plane_readv)(
                                    handle, "vplane", vplane,
                                    0, 0, ncomponents));
            }
            for (int j = 0; j < blocal; ++j) {
                for (int i = 0; i < alocal; ++i) {
                    for (int h = 0; h < ncomponents; ++h) {
                        const REAL expected =  (REAL) 5*((j + bstart) +  7)
                                                   +  2*((i + astart) +  3)
                                                   - h;
                        const REAL value = vplane[j*bstride + i*astride + h];
                        fct_chk_eq_dbl(value, expected);
                    }
                }
            }
            free(vplane);
            fct_req(0 == esio_file_close(handle));
        }
        FCT_TEST_END();

    }
    FCT_FIXTURE_SUITE_END();

    if (filetemplate) free(filetemplate);
}
FCT_END()
