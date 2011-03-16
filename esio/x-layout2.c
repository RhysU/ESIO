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

// TODO See bug #1229
// TODO See bug #1422

hid_t METHODNAME(hid_t plist_id, hid_t dset_id, QUALIFIER void *field,
                 int cglobal, int cstart, int clocal, int cstride,
                 int bglobal, int bstart, int blocal, int bstride,
                 int aglobal, int astart, int alocal, int astride,
                 hid_t type_id)
{
    (void) cglobal; /* Unused but present for API consistency */
    (void) aglobal; /* Unused but present for API consistency */

    /* Determine in-memory size of type_id */
    const size_t type_size = H5Tget_size(type_id);
    assert(type_size > 0);

    /* Temporaries used in hyperslab selection */
    hsize_t offset[2], stride[2], count[2];

    /* Establish (possibly strided) memspace details */
    count[0] = 1;
    count[1] = astride * alocal;
    const hid_t memspace = H5Screate_simple(2, count, NULL);
    assert(memspace > 0);
    offset[0] = 0;
    offset[1] = 0;
    stride[0] = bstride;
    stride[1] = astride;
    count[0]  = 1;
    count[1]  = alocal;
    if (H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
                            offset, stride, count, NULL) < 0) {
        H5Sclose(memspace);
        ESIO_ERROR("Selecting memory hyperslab failed", ESIO_EFAILED);
    }

    /* Establish filespace details */
    const hid_t filespace = H5Dget_space(dset_id);
    assert(filespace >= 0);

    for (int i = 0; i < clocal; ++i)
    {
        for (int j = 0; j < blocal; ++j)
        {
            /* Select hyperslab in the file */
            offset[0] = (j + bstart) + (i + cstart) * bglobal;
            offset[1] = astart;
            stride[0] = 1;
            stride[1] = 1;
            count[0]  = 1;
            count[1]  = alocal;
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
                                offset, stride, count, NULL);

            /* Compute memory offset to hyperslab's data */
            const size_t moffset = j*bstride + i*cstride;
#ifdef __INTEL_COMPILER
/* warning #1338: arithmetic on pointer to void or function type */
#pragma warning(push,disable:1338)
#endif
            /* Note use of type_size when adding to (void *) field */
            QUALIFIER void *p_field = field + (type_size * moffset);
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

            /* Transfer hyperslab to or from memory */
            const herr_t status = OPFUNC(dset_id, type_id, memspace,
                                         filespace, plist_id, p_field);
            if (status < 0) {
                H5Sclose(filespace);
                H5Sclose(memspace);
                ESIO_ERROR("Operation failed", ESIO_EFAILED);
            }
        }
    }

    /* Release temporary resources */
    H5Sclose(filespace);
    H5Sclose(memspace);

    return ESIO_SUCCESS;
}
