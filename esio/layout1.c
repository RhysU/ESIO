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

    /* Determine in-memory size of type_id */
    const size_t type_size = H5Tget_size(type_id);
    assert(type_size > 0);

    /* Establish contiguous filespace details */
    const hid_t filespace = H5Dget_space(dset_id);
    assert(filespace >= 0);

    /* Strategy changes depending on whether or not data is contiguous. */
    if (   (astride == 1)
        && (bstride == astride * alocal)
        && (cstride == bstride * blocal)) {

        /*
         * If contiguous, perform on single operation handling all data.
         */

        /* Establish contiguous memspace details */
        const hsize_t nelems = clocal * cstride;
        const hid_t memspace = H5Screate_simple(1, &nelems, NULL);
        assert(memspace > 0);

        /* Select appropriate hyperslab within the file */
        const hsize_t start[3] = { cstart, bstart, astart };
        const hsize_t count[3] = { clocal, blocal, alocal };
        if (H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
                                start, NULL, count, NULL) < 0) {
            H5Sclose(filespace);
            H5Sclose(memspace);
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

    } else {
        /*
         * Perform multiple regular hyperslab operations for strided memory.
         * See http://www.hdfgroup.org/HDF5/PHDF5/parallelhdf5hints.pdf
         * for the motivation.
         */

        /* Establish strided memspace details for a single pencil of data.   */
        /* We'll specify different in-memory offsets for this HDF5 concepts. */
        hsize_t count[3]  = { 1, 1, astride * alocal };
        const hid_t memspace = H5Screate_simple(3, count, NULL);
        assert(memspace > 0);
        hsize_t offset[3] = { 0, 0, 0                };
        hsize_t stride[3] = { 1, 1, astride          };
        count[0] = 1;
        count[1] = 1;
        count[2] = alocal;
        if (H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
                                offset, stride, count, NULL) < 0) {
            H5Sclose(filespace);
            H5Sclose(memspace);
            ESIO_ERROR("Selecting memory hyperslab failed", ESIO_EFAILED);
        }

        /* Loop over each pencil of data and perform the operation. */
        for (int i = 0; i < clocal; ++i) {
            for (int j = 0; j < blocal; ++j) {

                /* Select contiguous pencil within the file */
                offset[0] = i + cstart;
                offset[1] = j + bstart;
                offset[2] = astart;
                stride[0] = 1;
                stride[1] = 1;
                stride[2] = 1;
                count[0]  = 1;
                count[1]  = 1;
                count[2]  = alocal;
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
    }

    return ESIO_SUCCESS;
}
