/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of ESIO.
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
 * layout0.c: Implementations of field/vfield layout 0
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "error.h"
#include "layout.h"

hid_t esio_layout0_filespace_creator(int cglobal, int bglobal, int aglobal)
{
    const hsize_t dims[3] = { cglobal, bglobal, aglobal };
    return H5Screate_simple(3, dims, NULL);
}

// Reading and writing differ only by an HDF5 operation name and a const
// qualifier.  Define a macro that we'll use to implement both operations.
#define GEN_LAYOUT0_TRANSFER(METHODNAME, OPFUNC, QUALIFIER)                  \
hid_t METHODNAME(hid_t dset_id, QUALIFIER void *field,                       \
                 int cglobal, int cstart, int clocal, int cstride,           \
                 int bglobal, int bstart, int blocal, int bstride,           \
                 int aglobal, int astart, int alocal, int astride,           \
                 hid_t type_id)                                              \
{                                                                            \
    (void) cglobal; /* Unused but present for API consistency */             \
    (void) bglobal; /* Unused but present for API consistency */             \
    (void) aglobal; /* Unused but present for API consistency */             \
                                                                             \
    /* Create property list for collective operation */                      \
    const hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);                      \
    if (H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE) < 0) {              \
        H5Pclose(plist_id);                                                  \
        ESIO_ERROR(#METHODNAME ": setting IO transfer mode failed",          \
                   ESIO_EFAILED);                                            \
    }                                                                        \
                                                                             \
    /* Establish (possibly strided) memspace details */                      \
    const hsize_t nelems = clocal * cstride;                                 \
    const hid_t memspace = H5Screate_simple(1, &nelems, NULL);               \
    assert(memspace > 0);                                                    \
    for (int k = 0; k < clocal; ++k) {                                       \
        for (int j = 0; j < blocal; ++j) {                                   \
            const hsize_t start  = k*cstride + j*bstride;                    \
            const hsize_t stride = astride;                                  \
            const hsize_t count  = alocal;                                   \
            if (H5Sselect_hyperslab(memspace, H5S_SELECT_OR,                 \
                                    &start, &stride, &count, NULL) < 0) {    \
                H5Pclose(plist_id);                                          \
                H5Sclose(memspace);                                          \
                ESIO_ERROR(#METHODNAME ": selecting memory hyperslab failed",\
                           ESIO_EFAILED);                                    \
            }                                                                \
        }                                                                    \
    }                                                                        \
                                                                             \
    /* Establish (contiguous) filespace details */                           \
    const hid_t filespace = H5Dget_space(dset_id);                           \
    assert(filespace >= 0);                                                  \
    const hsize_t start[3] = { cstart, bstart, astart };                     \
    const hsize_t count[3] = { clocal, blocal, alocal };                     \
    if (H5Sselect_hyperslab(filespace, H5S_SELECT_SET,                       \
                            start, NULL, count, NULL) < 0) {                 \
        H5Pclose(plist_id);                                                  \
        H5Sclose(memspace);                                                  \
        H5Sclose(filespace);                                                 \
        ESIO_ERROR(#METHODNAME ": selecting file hyperslab failed",          \
                   ESIO_EFAILED);                                            \
    }                                                                        \
                                                                             \
    /* Transfer hyperslab to or from memory */                               \
    const herr_t status = OPFUNC(dset_id, type_id, memspace,                 \
                                 filespace, plist_id, field);                \
    if (status < 0) {                                                        \
        H5Sclose(filespace);                                                 \
        H5Sclose(memspace);                                                  \
        H5Pclose(plist_id);                                                  \
        ESIO_ERROR(#METHODNAME ": operation failed", ESIO_EFAILED);          \
    }                                                                        \
                                                                             \
    /* Release temporary resources */                                        \
    H5Sclose(filespace);                                                     \
    H5Sclose(memspace);                                                      \
    H5Pclose(plist_id);                                                      \
                                                                             \
    return ESIO_SUCCESS;                                                     \
}

// Routine to transfer data from const buffer to storage
GEN_LAYOUT0_TRANSFER(esio_layout0_field_writer, H5Dwrite, const)

// Routine to transfer data from storage to mutable buffer
GEN_LAYOUT0_TRANSFER(esio_layout0_field_reader, H5Dread, /* mutable */)
