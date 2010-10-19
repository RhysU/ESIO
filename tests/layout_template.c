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
 * layout_template.c: unit test template for ESIO field/vfield operations
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
        "--distribute",
        "-d",
        FCTCL_STORE_VALUE,
        "Choose which direction (one of c, b, or a) is distributed"
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
        "Sets the number of components used in vfield test"
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
    {
        "--auxstride-c",
        NULL,
        FCTCL_STORE_VALUE,
        "Auxiliary stride added to the c direction"
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
    int cglobal, cstart, clocal;
    int bglobal, bstart, blocal;
    int aglobal, astart, alocal;
    const char *distribute = fctcl_val2("--distribute", "a");
    if (fctstr_ieq(distribute,"c")) {
        // Data uniformly partitioned in the slow index
        clocal = partitioned_size;
        cglobal  = clocal * world_size;
        cstart = clocal * world_rank;

        aglobal  = alocal = unpartitioned_size;
        bglobal  = blocal = unpartitioned_size;
        astart = bstart = 0;
    } else if (fctstr_ieq(distribute,"b")) {
        // Data uniformly partitioned in the medium index
        blocal = partitioned_size;
        bglobal = blocal * world_size;
        bstart = blocal * world_rank;

        aglobal = alocal = unpartitioned_size;
        cglobal = clocal = unpartitioned_size;
        astart = cstart = 0;
    } else if (fctstr_ieq(distribute,"a")) {
        // Data uniformly partitioned in the fastest index
        alocal = partitioned_size;
        aglobal = alocal * world_size;
        astart = alocal * world_rank;

        bglobal = blocal = unpartitioned_size;
        cglobal = clocal = unpartitioned_size;
        bstart = cstart = 0;
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
    const int auxstride_c = (int) strtol(
        fctcl_val2("--auxstride-c","0"), (char **) NULL, 10);
    if (auxstride_c < 0) {
        fprintf(stderr, "\n--auxstride-c=%d < 0\n", auxstride_c);
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
    esio_state state;

    FCT_FIXTURE_SUITE_BGN(field_realvalued)
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
            state = esio_init(MPI_COMM_WORLD);
            assert(state);

            esio_layout_set(state, LAYOUT_TAG);
        }
        FCT_SETUP_END();

        FCT_TEARDOWN_BGN()
        {
            // Finalize ESIO state
            esio_finalize(state);

            // Clean up the unique file and filename
            if (world_rank == 0) {
                if (!preserve) unlink(filename);
            }
            if (filename) free(filename);
        }
        FCT_TEARDOWN_END();

        // Test scalar-valued fields, including overwrite details
        FCT_TEST_BGN(field)
        {
            fct_req(LAYOUT_TAG < esio_layout_count());
            fct_req(LAYOUT_TAG == esio_layout_get(state));

            TEST_REAL *field;

            // Compute stride information for each direction;
            // strides expressed in sizeof(TEST_REAL) which happens to be 1.
            const int astride  = 1                + auxstride_a * 1;
            const int bstride  = alocal * astride + auxstride_b * 1;
            const int cstride  = blocal * bstride + auxstride_c * 1;
            const size_t nelem = clocal * cstride;

            // Allocate storage for local portion of global field
            field = calloc(nelem, sizeof(TEST_REAL));
            fct_req(field);

            // Open file
            fct_req(0 == esio_file_create(state, filename, 1));

            // Write all zeros to disk
            fct_req(0 == TEST_ESIO_FIELD_WRITE(
                                state, "field", field,
                                cglobal, cstart, clocal, cstride,
                                bglobal, bstart, blocal, bstride,
                                aglobal, astart, alocal, astride));

            // Populate local field with test data
            for (int k = 0; k < clocal; ++k) {
                for (int j = 0; j < blocal; ++j) {
                    for (int i = 0; i < alocal; ++i) {
                        const TEST_REAL value = 11*((k + cstart) + 13)
                                              +  5*((j + bstart) +  7)
                                              +  2*((i + astart) +  3);
                        field[k*cstride + j*bstride + i*astride] = value;
                    }
                }
            }

            // Overwrite zeros on disk with test data
            if (auxstride_a || auxstride_b || auxstride_c) {
                // Noncontiguous; exercise non-default stride arguments
                fct_req(0 == TEST_ESIO_FIELD_WRITE(
                                    state, "field", field,
                                    cglobal, cstart, clocal, cstride,
                                    bglobal, bstart, blocal, bstride,
                                    aglobal, astart, alocal, astride));
            } else {
                // Contiguous; exercise default stride arguments
                fct_req(0 == TEST_ESIO_FIELD_WRITE(
                                    state, "field", field,
                                    cglobal, cstart, clocal, 0,
                                    bglobal, bstart, blocal, 0,
                                    aglobal, astart, alocal, 0));
            }

            { // Ensure the global size was written correctly
                int tmp_aglobal, tmp_bglobal, tmp_cglobal;
                fct_req(0 == esio_field_size(state, "field",
                                             &tmp_cglobal,
                                             &tmp_bglobal,
                                             &tmp_aglobal));
                fct_chk_eq_int(cglobal, tmp_cglobal);
                fct_chk_eq_int(bglobal, tmp_bglobal);
                fct_chk_eq_int(aglobal, tmp_aglobal);
            }

            // Close the file
            fct_req(0 == esio_file_close(state));

            // Free the field
            free(field);

            // Reopen the file using normal HDF5 APIs on root processor
            // Examine the contents (contiguously) to ensure it matches
            if (world_rank == 0) {
                field = calloc(cglobal * bglobal * aglobal,
                               sizeof(TEST_REAL));
                fct_req(field);
                const hid_t file_id
                    = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
                fct_req(0 <= H5LTread_dataset(file_id, "field",
                                              TEST_H5T, field));

                TEST_REAL *p_field = field;
                for (int k = 0; k < cglobal; ++k) {
                    for (int j = 0; j < bglobal; ++j) {
                        for (int i = 0; i < aglobal; ++i) {
                            const TEST_REAL expected
                                    = 2*(i+3)+5*(j+7)+11*(k+13);
                            const TEST_REAL value
                                    = *p_field++;
                            fct_chk_eq_dbl(value, expected);
                        }
                    }
                }

                fct_req(0 <= H5Fclose(file_id));
                free(field);
            }

            // Re-read the file in a distributed manner and verify contents
            // TODO Ensure non-referenced memory locations remain unmodified
            field = calloc(nelem, sizeof(TEST_REAL));
            fct_req(field);
            fct_req(0 == esio_file_open(state, filename, 0));
            if (auxstride_a || auxstride_b || auxstride_c) {
                // Noncontiguous; exercise non-default stride arguments
                fct_req(0 == TEST_ESIO_FIELD_READ(
                                    state, "field", field,
                                    cglobal, cstart, clocal, cstride,
                                    bglobal, bstart, blocal, bstride,
                                    aglobal, astart, alocal, astride));
            } else {
                // Contiguous; exercise default stride arguments
                fct_req(0 == TEST_ESIO_FIELD_READ(
                                    state, "field", field,
                                    cglobal, cstart, clocal, 0,
                                    bglobal, bstart, blocal, 0,
                                    aglobal, astart, alocal, 0));
            }
            for (int k = 0; k < clocal; ++k) {
                for (int j = 0; j < blocal; ++j) {
                    for (int i = 0; i < alocal; ++i) {
                        const TEST_REAL expected = 11*((k + cstart) + 13)
                                                 +  5*((j + bstart) +  7)
                                                 +  2*((i + astart) +  3);
                        const TEST_REAL value
                                = field[k*cstride + j*bstride + i*astride];
                        fct_chk_eq_dbl(value, expected);
                    }
                }
            }
            free(field);
            fct_req(0 == esio_file_close(state));
        }
        FCT_TEST_END();

        // Test vector-valued fields
        FCT_TEST_BGN(vfield)
        {
            fct_req(LAYOUT_TAG < esio_layout_count());
            fct_req(LAYOUT_TAG == esio_layout_get(state));

            TEST_REAL *vfield;

            // Compute stride information for each direction;
            // strides expressed in sizeof(TEST_REAL).
            const int astride  = ncomponents      + auxstride_a * ncomponents;
            const int bstride  = alocal * astride + auxstride_b * ncomponents;
            const int cstride  = blocal * bstride + auxstride_c * ncomponents;
            const size_t nelem = clocal * cstride;

            // Allocate storage for local portion of global vfield
            vfield = calloc(nelem, sizeof(TEST_REAL));
            fct_req(vfield);

            // Open file
            fct_req(0 == esio_file_create(state, filename, 1));

            // Populate local vfield with test data
            for (int k = 0; k < clocal; ++k) {
                for (int j = 0; j < blocal; ++j) {
                    for (int i = 0; i < alocal; ++i) {
                        for (int h = 0; h < ncomponents; ++h) {
                            vfield[k*cstride + j*bstride + i*astride + h]
                                    = 11*((k + cstart) + 13)
                                    +  5*((j + bstart) +  7)
                                    +  2*((i + astart) +  3)
                                    - h;
                        }
                    }
                }
            }

            // Write vfield to file
            if (auxstride_a || auxstride_b || auxstride_c) {
                // Noncontiguous; exercise non-default stride arguments
                fct_req(0 == TEST_ESIO_VFIELD_WRITE(
                                    state, "vfield", vfield,
                                    cglobal, cstart, clocal, cstride,
                                    bglobal, bstart, blocal, bstride,
                                    aglobal, astart, alocal, astride,
                                    ncomponents));
            } else {
                // Contiguous; exercise default stride arguments
                fct_req(0 == TEST_ESIO_VFIELD_WRITE(
                                    state, "vfield", vfield,
                                    cglobal, cstart, clocal, 0,
                                    bglobal, bstart, blocal, 0,
                                    aglobal, astart, alocal, 0,
                                    ncomponents));
            }

            { // Ensure the global size was written correctly
                int tmp_aglobal, tmp_bglobal, tmp_cglobal, tmp_ncomponents;
                fct_req(0 == esio_vfield_size(state, "vfield",
                                              &tmp_cglobal,
                                              &tmp_bglobal,
                                              &tmp_aglobal,
                                              &tmp_ncomponents));

                fct_chk_eq_int(cglobal,     tmp_cglobal);
                fct_chk_eq_int(bglobal,     tmp_bglobal);
                fct_chk_eq_int(aglobal,     tmp_aglobal);
                fct_chk_eq_int(ncomponents, tmp_ncomponents);
            }

            // Close the file
            fct_req(0 == esio_file_close(state));

            // Free the vfield
            free(vfield);

            // Reopen the file using normal HDF5 APIs on root processor
            // Examine the contents to ensure it matches
            if (world_rank == 0) {
                vfield = calloc(cglobal * bglobal * aglobal * ncomponents,
                                sizeof(TEST_REAL));
                fct_req(vfield);
                const hid_t file_id
                    = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
                hid_t type_id = TEST_H5T;
                if (ncomponents > 1) {
                    const hsize_t dims[1] = { ncomponents };
                    type_id = H5Tarray_create2(TEST_H5T, 1, dims);
                    fct_req(type_id >= 0);
                }
                fct_req(0 <= H5LTread_dataset(file_id, "vfield",
                                              type_id, vfield));
                if (ncomponents > 1) H5Tclose(type_id);

                TEST_REAL *p_field = vfield;
                for (int k = 0; k < cglobal; ++k) {
                    for (int j = 0; j < bglobal; ++j) {
                        for (int i = 0; i < aglobal; ++i) {
                            for (int h = 0; h < ncomponents; ++h) {
                                const TEST_REAL value = *p_field++;
                                fct_chk_eq_dbl(
                                    value,
                                    (TEST_REAL) 2*(i+3)+5*(j+7)+11*(k+13)-h);
                            }
                        }
                    }
                }

                fct_req(0 <= H5Fclose(file_id));
                free(vfield);
            }

            // Re-read the file in a distributed manner and verify contents
            // TODO Ensure non-referenced memory locations remain unmodified
            vfield = calloc(nelem, sizeof(TEST_REAL));
            fct_req(vfield);
            fct_req(0 == esio_file_open(state, filename, 0));
            if (auxstride_a || auxstride_b || auxstride_c) {
                // Noncontiguous; exercise non-default stride arguments
                fct_req(0 == TEST_ESIO_VFIELD_READ(
                                    state, "vfield", vfield,
                                    cglobal, cstart, clocal, cstride,
                                    bglobal, bstart, blocal, bstride,
                                    aglobal, astart, alocal, astride,
                                    ncomponents));
            } else {
                // Contiguous; exercise default stride arguments
                fct_req(0 == TEST_ESIO_VFIELD_READ(
                                    state, "vfield", vfield,
                                    cglobal, cstart, clocal, 0,
                                    bglobal, bstart, blocal, 0,
                                    aglobal, astart, alocal, 0,
                                    ncomponents));
            }
            for (int k = 0; k < clocal; ++k) {
                for (int j = 0; j < blocal; ++j) {
                    for (int i = 0; i < alocal; ++i) {
                        for (int h = 0; h < ncomponents; ++h) {
                            const TEST_REAL expected = 11*((k + cstart) + 13)
                                                     +  5*((j + bstart) +  7)
                                                     +  2*((i + astart) +  3)
                                                     - h;
                            const TEST_REAL value = vfield[
                                    k*cstride + j*bstride + i*astride + h];
                            fct_chk_eq_dbl(value, expected);
                        }
                    }
                }
            }
            free(vfield);
            fct_req(0 == esio_file_close(state));
        }
        FCT_TEST_END();

    }
    FCT_FIXTURE_SUITE_END();
}
FCT_END()
