//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// ExaScale IO library for turbulence simulation restart files
// http://github.com/RhysU/ESIO
//
// Copyright (C) 2010-2017, 2022, 2026 Rhys Ulerich
// Copyright (C) 2010-2017 The PECOS Development Team
//
// This file is part of ESIO.
//
// ESIO is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3.0 of the License, or
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "h5utils.h"

htri_t
esio_H5LTget_attribute_ndims_info(hid_t loc_id,
                                  const char *obj_name,
                                  const char *attr_name,
                                  int *rank,
                                  hsize_t *dims,
                                  H5T_class_t *type_class,
                                  size_t *type_size)
{
    // Blends H5LTget_attribute_ndims and H5LTget_attribute_info functions
    // Combining them allows reduced resource usage relative to separate calls
    // Attempt is made to preserve error codes returned by HDF5 API.

    htri_t result = -7;
    hid_t  obj_id  = H5I_INVALID_HID;
    hid_t  attr_id = H5I_INVALID_HID;
    hid_t  tid     = H5I_INVALID_HID;
    hid_t  sid     = H5I_INVALID_HID;
    H5T_class_t tmp_class;
    size_t tmp_size;
    int tmp_rank;

    /* Open the object */
    obj_id = H5Oopen(loc_id, obj_name, H5P_DEFAULT);
    if (obj_id == H5I_INVALID_HID) {
        result = 0;  /* Not found is not an error */
        goto bail;
    }

    /* Open the attribute. */
    attr_id = H5Aopen(obj_id, attr_name, H5P_DEFAULT);
    if (attr_id == H5I_INVALID_HID) {
        result = 0;  /* Not found is not an error */
        goto bail;
    }

    /* Get an identifier for the datatype. */
    tid = H5Aget_type(attr_id);
    if (tid == H5I_INVALID_HID) {
        result = -1;  /* Unexpected for valid (obj_id, attr_id) */
        goto bail;
    }

    /* Get the dataspace handle */
    sid = H5Aget_space(attr_id);
    if (sid == H5I_INVALID_HID) {
        result = -2;  /* Unexpected for valid (obj_id, attr_id) */
        goto bail;
    }

    /* Get the class without clobbering inputs until success */
    tmp_class = H5Tget_class(tid);
    if (tmp_class == H5T_NO_CLASS) {
        result = -3;  /* Unexpected for valid (obj_id, attr_id) */
        goto bail;
    }

    /* Get the size without clobbering inputs until success */
    tmp_size = H5Tget_size(tid);
    if (tmp_size == 0) {
        result = -4;  /* Unexpected for valid (obj_id, attr_id) */
        goto bail;
    }

    /* Get rank without clobbering inputs until success */
    tmp_rank = H5Sget_simple_extent_ndims(sid);
    if (tmp_rank < 0) {
        result = -5;  /* Unexpected for valid (obj_id, attr_id) */
        goto bail;
    }

    /* Get dimensions finally emitting outputs */
    if (H5Sget_simple_extent_dims(sid, dims, NULL) < 0) {
        result = -6;
        goto bail;
    }

    /* Now, overwrite other outputs */
    *type_class = tmp_class;
    *type_size = tmp_size;
    *rank = tmp_rank;
    result = 1;

bail:
    if (sid     != H5I_INVALID_HID) H5Sclose(sid);
    if (tid     != H5I_INVALID_HID) H5Tclose(tid);
    if (attr_id != H5I_INVALID_HID) H5Aclose(attr_id);
    if (obj_id  != H5I_INVALID_HID) H5Oclose(obj_id);
    return result;
}

struct qw_data {
    hid_t errnum;
    int   maj_pos;
    int   min_pos;
};

static herr_t
query_walker(unsigned n, const H5E_error2_t *err_desc, void *client_data)
{
    struct qw_data * const qw = (struct qw_data*) client_data;
    if (err_desc->maj_num == qw->errnum) qw->maj_pos = n;
    if (err_desc->min_num == qw->errnum) qw->min_pos = n;
    return 0;
}

htri_t
esio_H5Equery_stack(hid_t errnum)
{
    herr_t e = 0;

    /* Walk stack looking for errnum */
    struct qw_data qw = {errnum, -1, -1};
    if ((e = H5Ewalk2(H5E_DEFAULT, H5E_WALK_UPWARD, &query_walker, &qw)) < 0)
        return e;

    /* Return true iff we found the specified error code */
    return (qw.maj_pos == -1 && qw.min_pos == -1) ? 0 : 1;
}
