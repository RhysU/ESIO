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

#include <assert.h>
#include <limits.h>
#include <hdf5.h>
#include <hdf5_hl.h>

#include "error.h"
#include "esio.h"
#include "h5utils.h"
#include "metadata.h"
#include "version.h"

// TODO Have esio_XXX_metadata_read detect different sorts of errors

#define ESIO_FIELD_METADATA_SIZE (8)

static
int esio_hdf5metadata_read(hid_t loc_id,
                           const char *name,
                           const int rank,
                           int *global,
                           int *ncomponents);

int esio_type_ncomponents(hid_t type_id)
{
    hsize_t ncomponents = 0;

    // Look up number of scalar components contained in the given type
    switch (H5Tget_class(type_id)) {
        case H5T_ENUM:
        case H5T_FLOAT:
        case H5T_INTEGER:
        case H5T_OPAQUE:
            // Scalar types contain a single component for metadata purposes
            ncomponents = 1;
            break;
        case H5T_ARRAY:
            // One dimensional arrays have ncomponents equal to their size
            assert(H5Tget_array_ndims(type_id) == 1);
            H5Tget_array_dims2(type_id, &ncomponents);
            break;
        case H5T_COMPOUND:
            ESIO_ERROR_VAL("H5T_COMPOUND not supported", ESIO_ESANITY, -1);
        case H5T_REFERENCE:
            ESIO_ERROR_VAL("H5T_REFERENCE not supported", ESIO_ESANITY, -1);
        case H5T_STRING:
            ESIO_ERROR_VAL("H5T_STRING not supported", ESIO_ESANITY, -1);
        case H5T_VLEN:
            ESIO_ERROR_VAL("H5T_VLEN not supported", ESIO_ESANITY, -1);
        case H5T_TIME:
            ESIO_ERROR_VAL("H5T_TIME not supported", ESIO_ESANITY, -1);
        default:
            ESIO_ERROR_VAL("Unknown H5T_class_t value", ESIO_ESANITY, -1);
            break;
    }

    // Coerce to int and return
    assert(ncomponents > 0 && ncomponents <= INT_MAX);
    return (int) ncomponents;
}

hid_t esio_type_arrayify(hid_t type_id, int ncomponents)
{
    hid_t retval;

    if (ncomponents < 1) {
        ESIO_ERROR_VAL("ncomponents < 1", ESIO_EINVAL, -1);
    } else if (ncomponents == 1) {
        retval = H5Tcopy(type_id);
    } else {
        const hsize_t dims[1] = { ncomponents };
        retval = H5Tarray_create2(type_id, 1, dims);
    }

    if (retval < 0) {
        ESIO_ERROR_VAL("Error creating array of specified type",
                       ESIO_EFAILED, retval);
    }

    return retval;
}

int esio_field_metadata_write(hid_t loc_id, const char *name,
                              int layout_index,
                              int cglobal, int bglobal, int aglobal,
                              hid_t type_id)
{
    // layout0 is the default and requires writing no auxiliary metadata.
    // Requires companion logic in esio_field_metadata_read(...) below.
    if (layout_index == 0) return ESIO_SUCCESS;

    // Meant to be opaque but the cool kids will figure it out. :P
    const int ncomponents = esio_type_ncomponents(type_id);
    const int metadata[ESIO_FIELD_METADATA_SIZE] = {
        ESIO_MAJOR_VERSION,
        ESIO_MINOR_VERSION,
        ESIO_POINT_VERSION,
        layout_index,
        cglobal,
        bglobal,
        aglobal,
        ncomponents
    };
    const herr_t status = H5LTset_attribute_int(loc_id, name,
            "esio_field_metadata", metadata, ESIO_FIELD_METADATA_SIZE);
    return (status >= 0) ? ESIO_SUCCESS : ESIO_EFAILED;
}

int esio_field_metadata_read(hid_t loc_id, const char *name,
                             int *layout_index,
                             int *cglobal, int *bglobal, int *aglobal,
                             int *ncomponents)
{
    // This routine should not (generally) invoke any ESIO error handling
    // unless Bad Things (TM) happen.  It is sometimes used to query for the
    // existence of a field.

    // Local scratch space into which we read a field's metadata.
    // Employ a sentinel to balk if/when we accidentally blow out the buffer
    int metadata[ESIO_FIELD_METADATA_SIZE + 1];
    const int sentinel = INT_MIN + 999983;
    metadata[ESIO_FIELD_METADATA_SIZE] = sentinel;

    // Disable the error handler once and only once and perform all queries.
    // Dumb but required due to how disable/enable macros are implemented.
    DISABLE_HDF5_ERROR_HANDLER(one)
    const htri_t metadata_exists = H5Aexists_by_name(
            loc_id, name, "esio_field_metadata", H5P_DEFAULT);
    ENABLE_HDF5_ERROR_HANDLER(one)

    // If metadata did not exist use layout0 and employ esio_hdf5metadata_read.
    if (metadata_exists == 0) {
        assert(ESIO_FIELD_METADATA_SIZE >= 4);
        if (esio_hdf5metadata_read(loc_id, name, 3, metadata, metadata + 3)) {
            ESIO_ERROR("ESIO unable to read field (?) lacking esio_field_metadata",
                       ESIO_EFAILED); // Moderately Bad (TM)
        }
        if (layout_index) *layout_index = 0;
        if (cglobal)      *cglobal      = metadata[0];
        if (bglobal)      *bglobal      = metadata[1];
        if (aglobal)      *aglobal      = metadata[2];
        if (ncomponents)  *ncomponents  = metadata[3];
        return ESIO_SUCCESS;
    }

    // Metadata existed so read it into local scratch space
    DISABLE_HDF5_ERROR_HANDLER(two)
    const herr_t err = H5LTget_attribute_int(
            loc_id, name, "esio_field_metadata", metadata);
    ENABLE_HDF5_ERROR_HANDLER(two)

    // Check that our sentinel survived the read process (Very Bad (TM))
    if (metadata[ESIO_FIELD_METADATA_SIZE] != sentinel) {
        ESIO_ERROR("detected metadata buffer overflow", ESIO_ESANITY);
    }

    if (err < 0) {
        // On error, report "soft" failure to the caller
        return ESIO_EFAILED;
    } else {
        // On success, sanity check layout_index's value...
        if (metadata[3] < 0 || metadata[3] >= esio_field_layout_count()) {
            ESIO_ERROR("ESIO metadata contains unknown layout_index",
                       ESIO_ESANITY); // Very Bad (TM)
        }

        // ...and populate all requested, outgoing arguments.
        if (layout_index) *layout_index = metadata[3];
        if (cglobal)      *cglobal      = metadata[4];
        if (bglobal)      *bglobal      = metadata[5];
        if (aglobal)      *aglobal      = metadata[6];
        if (ncomponents)  *ncomponents  = metadata[7];

        return ESIO_SUCCESS;
    }
}

static
int esio_hdf5metadata_read(hid_t loc_id,
                           const char *name,
                           const int rank,
                           int *global,
                           int *ncomponents)
{
    // Extract metadata using HDF5's introspection utilities

    // Open dataset and invoke ESIO_ERROR for non-H5E_NOTFOUND errors.
    DISABLE_HDF5_ERROR_HANDLER(one)
    const hid_t dset_id = H5Dopen1(loc_id, name);
    ENABLE_HDF5_ERROR_HANDLER(one)
    if (dset_id == H5E_NOTFOUND) {
        return ESIO_NOTFOUND;
    } else if (dset_id < 0) {
        return ESIO_EINVAL;
    }

    // Obtain type information
    const hid_t type_id = H5Dget_type(dset_id);
    if (type_id < 0) {
        H5Dclose(dset_id);
        ESIO_ERROR("Unable to get datatype", ESIO_EFAILED);
    }

    // Obtain number of components from the type
    const int tmp_ncomponents = esio_type_ncomponents(type_id);
    if (tmp_ncomponents < 1) {
        H5Dclose(type_id);
        H5Dclose(dset_id);
        ESIO_ERROR("Unable to get datatype ncomponents", ESIO_EFAILED);
    }

    // Close type
    H5Tclose(type_id);

    // Open dataspace
    const hid_t space_id = H5Dget_space(dset_id);
    if (space_id < 0) {
        H5Dclose(dset_id);
        ESIO_ERROR("Unable to open dataspace", ESIO_ESANITY);
    }

    // Obtain the dataspace's extents
    hsize_t dims[H5S_MAX_RANK];
    if (rank != H5Sget_simple_extent_dims(space_id, dims, NULL)) {
        H5Sclose(space_id);
        H5Dclose(dset_id);
        ESIO_ERROR("Incorrect rank supplied for data", ESIO_EINVAL);
    }

    // Close resources
    H5Sclose(space_id);
    H5Dclose(dset_id);

    // Successfully retrieved all information; mutate arguments
    if (global) {
        for (int i = 0; i < rank; ++i)
            global[i] = dims[i];
    }
    if (ncomponents) *ncomponents = tmp_ncomponents;

    return ESIO_SUCCESS;
}

int esio_plane_metadata_write(hid_t loc_id, const char *name,
                              int bglobal, int aglobal,
                              hid_t type_id)
{
    (void) loc_id;  // Unused
    (void) name;    // Unused
    (void) bglobal; // Unused
    (void) aglobal; // Unused
    (void) type_id; // Unused
    return ESIO_SUCCESS; // NOP: Using esio_hdf5metadata_read to retrieve
}

int esio_plane_metadata_read(hid_t loc_id, const char *name,
                             int *bglobal, int *aglobal,
                             int *ncomponents)
{
    int global[2];
    const int status
        = esio_hdf5metadata_read(loc_id, name, 2, global, ncomponents);
    if (status == ESIO_SUCCESS) {
        if (bglobal) *bglobal = global[0];
        if (aglobal) *aglobal = global[1];
        // ncomponents modified in esio_hdf5metatadata_read
    }
    return status;
}

int esio_line_metadata_write(hid_t loc_id, const char *name,
                             int aglobal,
                             hid_t type_id)
{
    (void) loc_id;  // Unused
    (void) name;    // Unused
    (void) aglobal; // Unused
    (void) type_id; // Unused
    return ESIO_SUCCESS; // NOP: Using esio_hdf5metadata_read to retrieve
}

int esio_line_metadata_read(hid_t loc_id, const char *name,
                            int *aglobal,
                            int *ncomponents)
{
    return esio_hdf5metadata_read(loc_id, name, 1, aglobal, ncomponents);
}
