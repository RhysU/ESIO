//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// esio 0.1.3: ExaScale IO library for turbulence simulation restart files
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

    hid_t       attr_id;
    hid_t       tid;
    hid_t       sid;
    hid_t       obj_id;

    /* Open the object */
    if((obj_id = H5Oopen(loc_id, obj_name, H5P_DEFAULT)) < 0)
        return -1;

    /* Open the attribute. */
    if((attr_id = H5Aopen(obj_id, attr_name, H5P_DEFAULT)) < 0)
    {
        H5Oclose(obj_id);
        return -1;
    }

    /* Get an identifier for the datatype. */
    tid = H5Aget_type(attr_id);

    /* Get the class. */
    *type_class = H5Tget_class(tid);

    /* Get the size. */
    *type_size = H5Tget_size( tid );

    /* Get the dataspace handle */
    if ( (sid = H5Aget_space( attr_id )) < 0 )
        goto out;

    /* Get rank */
    if((*rank = H5Sget_simple_extent_ndims(sid)) < 0)
        goto out;

    /* Get dimensions */
    if ( H5Sget_simple_extent_dims( sid, dims, NULL) < 0 )
        goto out;

    /* Terminate access to the dataspace */
    if ( H5Sclose( sid ) < 0 )
        goto out;

    /* Release the datatype. */
    if ( H5Tclose( tid ) )
        goto out;

    /* End access to the attribute */
    if ( H5Aclose( attr_id ) )
        goto out;

    /* Close the object */
    if(H5Oclose(obj_id) < 0 )
        return -1;

    return 0;

out:
    H5Tclose(tid);
    H5Aclose(attr_id);
    H5Oclose(obj_id);
    return -1;
}
