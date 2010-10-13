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

hid_t esio_layout0_filespace_creator(int na, int nb, int nc)
{
    const hsize_t dims[2] = { nb * nc, na };
    return H5Screate_simple(2, dims, NULL);
}

hid_t esio_layout0_field_writer(hid_t dset_id, const void *field,
                                int na, int ast, int asz,
                                int nb, int bst, int bsz,
                                int nc, int cst, int csz,
                                hid_t type_id, size_t type_size)
{
    // TODO Error checking here

    (void) nc; // Unused but present for API consistency

    // Create property list for collective write
    const hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // Initialize one-time write details
    const hsize_t stride[2] = { 1, 1 };
    const hsize_t count[2]  = { 1, asz };
    const hid_t   memspace  = H5Screate_simple(2, count, NULL);
    const hid_t   filespace = H5Dget_space(dset_id);

    hsize_t offset[2];
    for (int i = 0; i < csz; ++i)
    {
        for (int j = 0; j < bsz; ++j)
        {
            // Select hyperslab in the file
            offset[0] = (j + bst) + (i + cst) * nb;
            offset[1] = ast;
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
                                offset, stride, count, NULL);

            // Compute in-memory offset to hyperslab's data
            // Note use of type_size when adding to (void *) field
            const size_t moffset = j*asz + i*asz*bsz;
#ifdef __INTEL_COMPILER
#pragma warning(push,disable:1338)
#endif
            const void *p_field = field + (type_size * moffset);
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

            // Write to hyperslab from memory
            const herr_t status = H5Dwrite(dset_id, type_id, memspace,
                                           filespace, plist_id, p_field);
            if (status < 0) {
                H5Sclose(filespace);
                H5Sclose(memspace);
                H5Pclose(plist_id);
                ESIO_ERROR("Write failed", ESIO_EFAILED);
            }
        }
    }

    // Release temporary resources
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);

    return ESIO_SUCCESS;
}

int esio_layout0_field_reader(hid_t dset_id, void *field,
                              int na, int ast, int asz,
                              int nb, int bst, int bsz,
                              int nc, int cst, int csz,
                              hid_t type_id, size_t type_size)
{
    return ESIO_ESANITY; // FIXME Implement
}
