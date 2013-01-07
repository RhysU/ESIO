//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// ESIO 0.1.8: ExaScale IO library for turbulence simulation restart files
// http://red.ices.utexas.edu/projects/esio/
//
// Copyright (C) 2010, 2011, 2012, 2013 The PECOS Development Team
//
// This file is part of ESIO.
//
// ESIO is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 2.1 of the License, or
// (at your option) any later version.
//
// ESIO is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with ESIO.  If not, see <http://www.gnu.org/licenses/>.
//
//-----------------------------------------------------------------------el-
// $Id$

// Designed to be #included from layout.c
#if !defined(METHODNAME) || !defined(OPFUNC) || !defined(QUALIFIER)
#error "One of METHODNAME, OPFUNC, or QUALIFIER not defined"
#endif

hid_t METHODNAME(hid_t plist_id, hid_t dset_id, QUALIFIER void *field,
                 int cglobal, int cstart, int clocal, int cstride,
                 int bglobal, int bstart, int blocal, int bstride,
                 int aglobal, int astart, int alocal, int astride,
                 hid_t type_id)
{
    (void) cglobal; /* Unused but present for API consistency */
    (void) bglobal; /* Unused but present for API consistency */
    (void) aglobal; /* Unused but present for API consistency */

    /* Establish (possibly strided) memspace details */
    const hsize_t nelems = clocal * cstride;
    const hsize_t lies   = 1;
    const hid_t memspace = H5Screate_simple(
            1, clocal * blocal * alocal > 0 ? &nelems : &lies, NULL);
    assert(memspace > 0);
    if (clocal * blocal * alocal == 0) {
        H5Sselect_none(memspace);
    } else if (   astride != 1
               || bstride != astride * alocal
               || cstride != bstride * blocal) {
        /* Strided memspace; additional hyperslab selection necessary */
        if (H5Sselect_none(memspace) < 0) {
            H5Sclose(memspace);
            ESIO_ERROR("Resetting memory hyperslab failed", ESIO_EFAILED);
        }
        for (int k = 0; k < clocal; ++k) {
            for (int j = 0; j < blocal; ++j) {
                const hsize_t start  = k*cstride + j*bstride;
                const hsize_t stride = astride;
                const hsize_t count  = alocal;
                if (H5Sselect_hyperslab(memspace, H5S_SELECT_OR,
                                        &start, &stride, &count, NULL) < 0) {
                    H5Sclose(memspace);
                    ESIO_ERROR("Selecting memory hyperslab failed",
                               ESIO_EFAILED);
                }
            }
        }
    }

    /* Establish (contiguous) filespace details */
    const hid_t filespace = H5Dget_space(dset_id);
    assert(filespace >= 0);
    const hsize_t start[3] = { cstart, bstart, astart };
    const hsize_t count[3] = { clocal, blocal, alocal };
    if (clocal * blocal * alocal == 0) {
        H5Sselect_none(filespace);
    } else if (H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
                                   start, NULL, count, NULL) < 0) {
        H5Sclose(memspace);
        H5Sclose(filespace);
        ESIO_ERROR("Selecting file hyperslab failed", ESIO_EFAILED);
    }

    /* Transfer hyperslab to or from memory */
    const herr_t status = OPFUNC(dset_id, type_id, memspace,
                                 filespace, plist_id, field);
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
