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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <hdf5.h>
#include "h5utils.h"

herr_t
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

    hid_t  attr_id;
    hid_t  tid;
    hid_t  sid;
    hid_t  obj_id;
    herr_t e = 0;

    /* Open the object */
    if ((obj_id = H5Oopen(loc_id, obj_name, H5P_DEFAULT)) < 0)
        return obj_id;

    /* Open the attribute. */
    if ((e = attr_id = H5Aopen(obj_id, attr_name, H5P_DEFAULT)) < 0) {
        /* H5Aopen gives no information about existence as of HDF5 1.8.7 */
        if (H5Aexists(obj_id, attr_name) == 0)
            e = H5E_NOTFOUND;

        H5Oclose(obj_id);
        return e;
    }

    /* Get an identifier for the datatype. */
    tid = H5Aget_type(attr_id);

    /* Get the class. */
    *type_class = H5Tget_class(tid);

    /* Get the size. */
    *type_size = H5Tget_size(tid);

    /* Get the dataspace handle */
    if ((e = sid = H5Aget_space(attr_id)) < 0)
        goto bail;

    /* Get rank */
    if ((e = *rank = H5Sget_simple_extent_ndims(sid)) < 0)
        goto bail;

    /* Get dimensions */
    if ((e = H5Sget_simple_extent_dims(sid, dims, NULL)) < 0)
        goto bail;

    /* Terminate access to the dataspace */
    if ((e = H5Sclose(sid)) < 0)
        goto bail;

    /* Release the datatype. */
    if ((e = H5Tclose(tid)))
        goto bail;

    /* End access to the attribute */
    if ((e = H5Aclose(attr_id)) < 0)
        goto bail;

    /* Close the object */
    if ((e = H5Oclose(obj_id)) < 0)
        return e;

    return 0;

bail:
    H5Tclose(tid);
    H5Aclose(attr_id);
    H5Oclose(obj_id);
    return e ? e : -1;
}
