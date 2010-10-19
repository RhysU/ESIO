//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// ESIO - ExaScale IO library
//
// ESIO is a restart file library for HPC exascale IO geared especially for
// turbulence applications.
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <limits.h>
#include <hdf5.h>
#include <hdf5_hl.h>

#include "error.h"
#include "esio.h"
#include "version.h"
#include "metadata.h"

#define ESIO_FIELD_METADATA_SIZE (8)
#define ESIO_PLANE_METADATA_SIZE (6)
#define ESIO_LINE_METADATA_SIZE  (5)
#define ESIO_POINT_METADATA_SIZE (4)

// Macro to save and disable current HDF5 error handler
#define DISABLE_HDF5_ERROR_HANDLER                               \
    H5E_auto2_t hdf5_handler;                                    \
    void *hdf5_client_data;                                      \
    H5Eget_auto2(H5E_DEFAULT, &hdf5_handler, &hdf5_client_data); \
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

// Macro to restore previously saved HDF5 error handler
#define ENABLE_HDF5_ERROR_HANDLER                               \
    H5Eset_auto2(H5E_DEFAULT, hdf5_handler, hdf5_client_data);

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

herr_t esio_field_metadata_write(hid_t loc_id, const char *name,
                                 int layout_tag,
                                 int cglobal, int bglobal, int aglobal,
                                 hid_t type_id)
{
    // Meant to be opaque but the cool kids will figure it out. :P
    const int ncomponents = esio_type_ncomponents(type_id);
    const int metadata[ESIO_FIELD_METADATA_SIZE] = {
        ESIO_MAJOR_VERSION,
        ESIO_MINOR_VERSION,
        ESIO_POINT_VERSION,
        layout_tag,
        cglobal,
        bglobal,
        aglobal,
        ncomponents
    };
    return H5LTset_attribute_int(loc_id, name, "esio_field_metadata",
                                 metadata, ESIO_FIELD_METADATA_SIZE);
}

herr_t esio_field_metadata_read(hid_t loc_id, const char *name,
                                int *layout_tag,
                                int *cglobal, int *bglobal, int *aglobal,
                                int *ncomponents)
{
    // This routine should not invoke any ESIO error handling--
    // It is sometimes used to query for the existence of a field.

    // Obtain current HDF5 error handler
    H5E_auto2_t hdf5_handler;
    void *hdf5_client_data;
    H5Eget_auto2(H5E_DEFAULT, &hdf5_handler, &hdf5_client_data);

    // Disable HDF5 error handler during metadata read
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    // Local scratch space into which we read the metadata
    // Employ a sentinel to balk if/when we accidentally blow out the buffer
    int metadata[ESIO_FIELD_METADATA_SIZE + 1];
    const int sentinel = INT_MIN + 999983;
    metadata[ESIO_FIELD_METADATA_SIZE] = sentinel;

    // Read the metadata into the buffer
    const herr_t err = H5LTget_attribute_int(
            loc_id, name, "esio_field_metadata", metadata);

    // Re-enable the HDF5 error handler
    H5Eset_auto2(H5E_DEFAULT, hdf5_handler, hdf5_client_data);

    // Check that our sentinel survived the read process
    if (metadata[ESIO_FIELD_METADATA_SIZE] != sentinel) {
        ESIO_ERROR_VAL("detected metadata buffer overflow", ESIO_ESANITY, err);
    }

    // On success...
    if (err >= 0) {
        // ... populate all requested, outgoing arguments
        if (layout_tag)  *layout_tag  = metadata[3];
        if (cglobal)     *cglobal     = metadata[4];
        if (bglobal)     *bglobal     = metadata[5];
        if (aglobal)     *aglobal     = metadata[6];
        if (ncomponents) *ncomponents = metadata[7];

        // ... sanity check layout_tag's value
        if (metadata[3] < 0 || metadata[3] >= esio_layout_count()) {
            ESIO_ERROR_VAL("ESIO metadata contains unknown layout_tag",
                           ESIO_ESANITY, err);
        }
    }

    // Return the H5LTget_attribute_int error code
    return err;
}

herr_t esio_plane_metadata_write(hid_t loc_id, const char *name,
                                 int bglobal, int aglobal,
                                 hid_t type_id)
{
    const int ncomponents = esio_type_ncomponents(type_id);
    const int metadata[ESIO_PLANE_METADATA_SIZE] = {
        ESIO_MAJOR_VERSION,
        ESIO_MINOR_VERSION,
        ESIO_POINT_VERSION,
        bglobal,
        aglobal,
        ncomponents
    };
    return H5LTset_attribute_int(loc_id, name, "esio_plane_metadata",
                                 metadata, ESIO_PLANE_METADATA_SIZE);
}

herr_t esio_plane_metadata_read(hid_t loc_id, const char *name,
                                int *bglobal, int *aglobal,
                                int *ncomponents)
{
    // This routine should not invoke any ESIO error handling--
    // It is sometimes used to query for the existence of a plane.

    // Local scratch space into which we read the metadata
    // Employ a sentinel to balk if/when we accidentally blow out the buffer
    int metadata[ESIO_PLANE_METADATA_SIZE + 1];
    const int sentinel = INT_MIN + 999983;
    metadata[ESIO_PLANE_METADATA_SIZE] = sentinel;


    // Read the metadata into the buffer
    DISABLE_HDF5_ERROR_HANDLER
    const herr_t err = H5LTget_attribute_int(
            loc_id, name, "esio_plane_metadata", metadata);
    ENABLE_HDF5_ERROR_HANDLER

    // Check that our sentinel survived the read process
    if (metadata[ESIO_PLANE_METADATA_SIZE] != sentinel) {
        ESIO_ERROR_VAL("detected metadata buffer overflow", ESIO_ESANITY, err);
    }

    // On success populate all requested, outgoing arguments
    if (err >= 0) {
        if (bglobal)     *bglobal     = metadata[3];
        if (aglobal)     *aglobal     = metadata[4];
        if (ncomponents) *ncomponents = metadata[5];
    }

    // Return the H5LTget_attribute_int error code
    return err;
}

herr_t esio_line_metadata_write(hid_t loc_id, const char *name,
                                int aglobal,
                                hid_t type_id)
{
    const int ncomponents = esio_type_ncomponents(type_id);
    const int metadata[ESIO_LINE_METADATA_SIZE] = {
        ESIO_MAJOR_VERSION,
        ESIO_MINOR_VERSION,
        ESIO_POINT_VERSION,
        aglobal,
        ncomponents
    };
    return H5LTset_attribute_int(loc_id, name, "esio_line_metadata",
                                 metadata, ESIO_LINE_METADATA_SIZE);
}

herr_t esio_line_metadata_read(hid_t loc_id, const char *name,
                               int *aglobal,
                               int *ncomponents)
{
    // This routine should not invoke any ESIO error handling--
    // It is sometimes used to query for the existence of a line.

    // Local scratch space into which we read the metadata
    // Employ a sentinel to balk if/when we accidentally blow out the buffer
    int metadata[ESIO_LINE_METADATA_SIZE + 1];
    const int sentinel = INT_MIN + 999983;
    metadata[ESIO_LINE_METADATA_SIZE] = sentinel;

    // Read the metadata into the buffer
    DISABLE_HDF5_ERROR_HANDLER
    const herr_t err = H5LTget_attribute_int(
            loc_id, name, "esio_line_metadata", metadata);
    ENABLE_HDF5_ERROR_HANDLER

    // Check that our sentinel survived the read process
    if (metadata[ESIO_LINE_METADATA_SIZE] != sentinel) {
        ESIO_ERROR_VAL("detected metadata buffer overflow", ESIO_ESANITY, err);
    }

    // On success populate all requested, outgoing arguments
    if (err >= 0) {
        if (aglobal)     *aglobal     = metadata[3];
        if (ncomponents) *ncomponents = metadata[4];
    }

    // Return the H5LTget_attribute_int error code
    return err;
}

herr_t esio_point_metadata_write(hid_t loc_id, const char *name,
                                 hid_t type_id)
{
    const int ncomponents = esio_type_ncomponents(type_id);
    const int metadata[ESIO_POINT_METADATA_SIZE] = {
        ESIO_MAJOR_VERSION,
        ESIO_MINOR_VERSION,
        ESIO_POINT_VERSION,
        ncomponents
    };
    return H5LTset_attribute_int(loc_id, name, "esio_point_metadata",
                                 metadata, ESIO_POINT_METADATA_SIZE);
}

herr_t esio_point_metadata_read(hid_t loc_id, const char *name,
                                int *ncomponents)
{
    // This routine should not invoke any ESIO error handling--
    // It is sometimes used to query for the existence of a point.

    // Local scratch space into which we read the metadata
    // Employ a sentinel to balk if/when we accidentally blow out the buffer
    int metadata[ESIO_POINT_METADATA_SIZE + 1];
    const int sentinel = INT_MIN + 999983;
    metadata[ESIO_POINT_METADATA_SIZE] = sentinel;

    // Read the metadata into the buffer
    DISABLE_HDF5_ERROR_HANDLER
    const herr_t err = H5LTget_attribute_int(
            loc_id, name, "esio_point_metadata", metadata);
    ENABLE_HDF5_ERROR_HANDLER

    // Check that our sentinel survived the read process
    if (metadata[ESIO_POINT_METADATA_SIZE] != sentinel) {
        ESIO_ERROR_VAL("detected metadata buffer overflow", ESIO_ESANITY, err);
    }

    // On success populate all requested, outgoing arguments
    if (err >= 0) {
        if (ncomponents) *ncomponents = metadata[3];
    }

    // Return the H5LTget_attribute_int error code
    return err;
}
