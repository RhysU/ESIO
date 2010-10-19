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
 * layout.c: Implementations of various field layout options
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

// **************************************************************
// LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1
// **************************************************************

hid_t esio_layout0_filespace_creator(int cglobal, int bglobal, int aglobal)
{
    const hsize_t dims[2] = { cglobal * bglobal, aglobal };
    return H5Screate_simple(2, dims, NULL);
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
    /* Determine in-memory size of type_id */                                \
    const size_t type_size = H5Tget_size(type_id);                           \
    assert(type_size > 0);                                                   \
                                                                             \
    /* Temporaries used in hyperslab selection */                            \
    hsize_t offset[2], stride[2], count[2];                                  \
                                                                             \
    /* Establish (possibly strided) memspace details */                      \
    count[0] = 1;                                                            \
    count[1] = astride * alocal;                                             \
    const hid_t memspace = H5Screate_simple(2, count, NULL);                 \
    assert(memspace > 0);                                                    \
    offset[0] = 0;                                                           \
    offset[1] = 0;                                                           \
    stride[0] = bstride;                                                     \
    stride[1] = astride;                                                     \
    count[0]  = 1;                                                           \
    count[1]  = alocal;                                                      \
    if (H5Sselect_hyperslab(memspace, H5S_SELECT_SET,                        \
                            offset, stride, count, NULL) < 0) {              \
        H5Pclose(plist_id);                                                  \
        H5Sclose(memspace);                                                  \
        ESIO_ERROR(#METHODNAME ": selecting memory hyperslab failed",        \
                   ESIO_EFAILED);                                            \
    }                                                                        \
                                                                             \
    /* Establish filespace details */                                        \
    const hid_t filespace = H5Dget_space(dset_id);                           \
    assert(filespace >= 0);                                                  \
                                                                             \
    for (int i = 0; i < clocal; ++i)                                         \
    {                                                                        \
        for (int j = 0; j < blocal; ++j)                                     \
        {                                                                    \
            /* Select hyperslab in the file */                               \
            offset[0] = (j + bstart) + (i + cstart) * bglobal;               \
            offset[1] = astart;                                              \
            stride[0] = 1;                                                   \
            stride[1] = 1;                                                   \
            count[0]  = 1;                                                   \
            count[1]  = alocal;                                              \
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET,                   \
                                offset, stride, count, NULL);                \
                                                                             \
            /* Compute memory offset to hyperslab's data */                  \
            /* Note use of type_size when adding to (void *) field */        \
            const size_t moffset = j*bstride + i*cstride;                    \
            QUALIFIER void *p_field = field + (type_size * moffset);         \
                                                                             \
            /* Transfer hyperslab to or from memory */                       \
            const herr_t status = OPFUNC(dset_id, type_id, memspace,         \
                                         filespace, plist_id, p_field);      \
            if (status < 0) {                                                \
                H5Sclose(filespace);                                         \
                H5Sclose(memspace);                                          \
                H5Pclose(plist_id);                                          \
                ESIO_ERROR(#METHODNAME ": operation failed", ESIO_EFAILED);  \
            }                                                                \
        }                                                                    \
    }                                                                        \
                                                                             \
    /* Release temporary resources */                                        \
    H5Sclose(filespace);                                                     \
    H5Sclose(memspace);                                                      \
    H5Pclose(plist_id);                                                      \
                                                                             \
    return ESIO_SUCCESS;                                                     \
}


#ifdef __INTEL_COMPILER
/* warning #1338: arithmetic on pointer to void or function type */
#pragma warning(push,disable:1338)
#endif

// Routine to transfer data from const buffer to storage
GEN_LAYOUT0_TRANSFER(esio_layout0_field_writer, H5Dwrite, const)

// Routine to transfer data from storage to mutable buffer
GEN_LAYOUT0_TRANSFER(esio_layout0_field_reader, H5Dread, /* mutable */)

#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif
