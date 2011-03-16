//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// esio 0.1.6: ExaScale IO library for turbulence simulation restart files
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

// Designed to be #included from layout.c
#if !defined(METHODNAME) || !defined(OPFUNC) || !defined(QUALIFIER)
#error "One of METHODNAME, OPFUNC, or QUALIFIER not defined"
#endif

hid_t METHODNAME(hid_t plist_id, hid_t dset_id, QUALIFIER void *plane,
                 int bglobal, int bstart, int blocal, int bstride,
                 int aglobal, int astart, int alocal, int astride,
                 hid_t type_id)
{
    (void) bglobal; /* Unused but present for API consistency */
    (void) aglobal; /* Unused but present for API consistency */

    /* Establish (possibly strided) memspace details */
    const hsize_t nelems = blocal * bstride;
    const hsize_t lies   = 1;
    const hid_t memspace
        = H5Screate_simple(1, blocal * alocal ? &nelems : &lies, NULL);
    assert(memspace > 0);
    if (blocal * alocal == 0) {
        H5Sselect_none(memspace);
    } else if (astride != 1 || bstride != astride * alocal) {
        /* Strided memspace; additional hyperslab selection necessary */
        if (H5Sselect_none(memspace) < 0) {
            H5Sclose(memspace);
            ESIO_ERROR("Resetting memory hyperslab failed", ESIO_EFAILED);
        }
        for (int j = 0; j < blocal; ++j) {
            const hsize_t start  = j*bstride;
            const hsize_t stride = astride;
            const hsize_t count  = alocal;
            if (H5Sselect_hyperslab(memspace, H5S_SELECT_OR,
                                    &start, &stride, &count, NULL) < 0) {
                H5Sclose(memspace);
                ESIO_ERROR("Selecting memory hyperslab failed", ESIO_EFAILED);
            }
        }
    }

    /* Establish (contiguous) filespace details */
    const hid_t filespace = H5Dget_space(dset_id);
    assert(filespace >= 0);
    const hsize_t start[2] = { bstart, astart };
    const hsize_t count[2] = { blocal, alocal };
    if (blocal * alocal == 0) {
        H5Sselect_none(filespace);
    } else if (H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
                                   start, NULL, count, NULL) < 0) {
        H5Sclose(memspace);
        H5Sclose(filespace);
        ESIO_ERROR("Selecting file hyperslab failed", ESIO_EFAILED);
    }

    /* Transfer hyperslab to or from memory */
    const herr_t status = OPFUNC(dset_id, type_id, memspace,
                                 filespace, plist_id, plane);
    if (status < 0) {
        H5Sclose(filespace);
        H5Sclose(memspace);
        ESIO_ERROR("Operation failed", ESIO_EFAILED);
    }

    /* Release temporary resources */
    H5Sclose(filespace);
    H5Sclose(memspace);

    return ESIO_SUCCESS;
}
