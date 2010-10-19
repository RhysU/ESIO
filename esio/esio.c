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
 * esio.c: Primary implementation details
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <mpi.h>

#include "error.h"
#include "esio.h"
#include "layout.h"
#include "version.h"

//*********************************************************************
// INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL
//*********************************************************************

static
MPI_Comm esio_MPI_Comm_dup_with_name(MPI_Comm comm);

static
int esio_type_ncomponents(hid_t type_id);

static
hid_t esio_type_arrayify(hid_t type_id, int ncomponents);

static
herr_t esio_field_metadata_write(hid_t loc_id, const char *name,
                                 int layout_tag,
                                 int cglobal, int bglobal, int aglobal,
                                 int ncomponents);

static
herr_t esio_field_metadata_read(hid_t loc_id, const char *name,
                                int *layout_tag,
                                int *cglobal, int *bglobal, int *aglobal,
                                int *ncomponents);

// See esio_field_metadata_write and esio_field_metadata_read
#define ESIO_FIELD_METADATA_SIZE (8)

static
hid_t esio_field_create(const esio_state s,
                        int cglobal, int bglobal, int aglobal,
                        const char* name, hid_t type_id);

static
int esio_field_close(hid_t dataset_id);

static
int esio_field_write_internal(const esio_state s,
                              const char* name,
                              const void *field,
                              int cglobal, int cstart, int clocal, int cstride,
                              int bglobal, int bstart, int blocal, int bstride,
                              int aglobal, int astart, int alocal, int astride,
                              hid_t type_id);

static
int esio_field_read_internal(const esio_state s,
                             const char* name,
                             void *field,
                             int cglobal, int cstart, int clocal, int cstride,
                             int bglobal, int bstart, int blocal, int bstride,
                             int aglobal, int astart, int alocal, int astride,
                             hid_t type_id);

//*********************************************************************
// INTERNAL TYPES INTERNAL TYPES INTERNAL TYPES INTERNAL TYPES INTERNAL
//*********************************************************************

struct esio_state_s {
    MPI_Comm  comm;       //< Communicator used for collective calls
    int       comm_rank;  //< Process rank within in MPI communicator
    int       comm_size;  //< Number of ranks within MPI communicator
    MPI_Info  info;       //< Info object used for collective calls
    hid_t     file_id;    //< Active HDF file identifier
    int       layout_tag; //< Active field layout_tag within HDF5 file
};

//***************************************************************************
// IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION
//***************************************************************************

// Lookup table of all the different layout types we understand
static const struct {
    int                      tag;
    esio_filespace_creator_t filespace_creator;
    esio_field_writer_t      field_writer;
    esio_field_reader_t      field_reader;
} esio_layout[] = {
    {
        0,
        &esio_layout0_filespace_creator,
        &esio_layout0_field_writer,
        &esio_layout0_field_reader
    }
};
static const int esio_nlayout = sizeof(esio_layout)/sizeof(esio_layout[0]);


static
MPI_Comm
esio_MPI_Comm_dup_with_name(MPI_Comm comm)
{
    if (comm == MPI_COMM_NULL) return MPI_COMM_NULL;

#ifdef MPI_MAX_OBJECT_NAME
    char buffer[MPI_MAX_OBJECT_NAME] = "";
#else
    char buffer[MPI_MAX_NAME_STRING] = "";
#endif
    int resultlen = 0;
    const int get_name_error = MPI_Comm_get_name(comm, buffer, &resultlen);
    if (get_name_error) {
        ESIO_MPICHKR(get_name_error /* MPI_Comm_get_name */);
        return MPI_COMM_NULL;
    }

    MPI_Comm retval = MPI_COMM_NULL;

    const int dup_error = MPI_Comm_dup(comm, &retval);
    if (dup_error) {
        ESIO_MPICHKR(dup_error /* MPI_Comm_dup */);
        return MPI_COMM_NULL;
    }

    if (resultlen > 0) {
        const int set_name_error = MPI_Comm_set_name(retval, buffer);
        if (set_name_error) {
            ESIO_MPICHKR(set_name_error /* MPI_Comm_set_name */);
            ESIO_MPICHKR(MPI_Comm_free(&retval));
            return MPI_COMM_NULL;
        }
    }

    return retval;
}

esio_state
esio_init(MPI_Comm comm)
{
    // Sanity check incoming arguments
    if (comm == MPI_COMM_NULL) {
        ESIO_ERROR_NULL("comm == MPI_COMM_NULL required", ESIO_EINVAL);
    }

    // Get number of processors and the local rank within the communicator
    int comm_size;
    ESIO_MPICHKN(MPI_Comm_size(comm, &comm_size));
    int comm_rank;
    ESIO_MPICHKN(MPI_Comm_rank(comm, &comm_rank));

    // Initialize an MPI Info instance
    MPI_Info info;
    ESIO_MPICHKN(MPI_Info_create(&info));

    // Create and initialize ESIO's opaque state struct
    esio_state s = calloc(1, sizeof(struct esio_state_s));
    if (s == NULL) {
        ESIO_ERROR_NULL("failed to allocate space for state", ESIO_ENOMEM);
    }
    s->comm        = esio_MPI_Comm_dup_with_name(comm);
    s->comm_rank   = comm_rank;
    s->comm_size   = comm_size;
    s->info        = info;
    s->file_id     = -1;
    s->layout_tag  = 0;

    if (s->comm == MPI_COMM_NULL) {
        esio_finalize(s);
        ESIO_ERROR_NULL("Detected MPI_COMM_NULL in s->comm", ESIO_ESANITY);
    }

    return s;
}

#ifdef __INTEL_COMPILER
// remark #1418: external function definition with no prior declaration
#pragma warning(push,disable:1418)
#endif
esio_state
esio_init_fortran(MPI_Fint fcomm)
{
    // Converting MPI communicators from Fortran to C requires MPI_Comm_f2c
    // See section 16.3.4 of the MPI 2.2 Standard for details
    return esio_init(MPI_Comm_f2c(fcomm));
}
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

int
esio_finalize(esio_state s)
{
    if (s) {
        if (s->file_id != -1) {
            esio_file_close(s); // Force file closure
        }
        if (s->comm != MPI_COMM_NULL) {
            ESIO_MPICHKR(MPI_Comm_free(&s->comm));
            s->comm = MPI_COMM_NULL;
        }
        if (s->info != MPI_INFO_NULL) {
            ESIO_MPICHKR(MPI_Info_free(&s->info));
            s->info = MPI_INFO_NULL;
        }
        free(s);
    }

    return ESIO_SUCCESS;
}

int
esio_layout_get(const esio_state s)
{
    if (s == NULL) {
        ESIO_ERROR("s == NULL", ESIO_EINVAL);
    }

    return s->layout_tag;
}

int
esio_layout_set(esio_state s, int layout_tag)
{
    if (s == NULL) {
        ESIO_ERROR("s == NULL", ESIO_EINVAL);
    }
    if (layout_tag < 0) {
        ESIO_ERROR("layout_tag < 0", ESIO_EINVAL);
    }
    if (layout_tag >= esio_nlayout) {
        ESIO_ERROR("layout_tag >= esio_nlayout", ESIO_EINVAL);
    }

    s->layout_tag = layout_tag;

    return ESIO_SUCCESS;
}

int
esio_file_create(esio_state s, const char *file, int overwrite)
{
    // Sanity check incoming arguments
    if (s == NULL) {
        ESIO_ERROR("s == NULL", ESIO_EINVAL);
    }
    if (s->file_id != -1) {
        ESIO_ERROR("Cannot create file because previous file not closed",
                   ESIO_EINVAL);
    }
    if (file == NULL) {
        ESIO_ERROR("file == NULL", ESIO_EINVAL);
    }

    // Initialize file creation property list identifier
    const hid_t fcpl_id = H5P_DEFAULT;

    // Initialize file access list property identifier
    const hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    if (fapl_id == -1) {
        ESIO_ERROR("Unable to create fapl_id", ESIO_ESANITY);
    }
    // Set collective details
    if (H5Pset_fapl_mpio(fapl_id, s->comm, s->info)) {
        H5Pclose(fapl_id);
        ESIO_ERROR("Unable to store MPI details in fapl_id", ESIO_ESANITY);
    }

    // Collectively create the file
    hid_t file_id;
    if (overwrite) {
        // Create new file or truncate existing file
        file_id = H5Fcreate(file, H5F_ACC_TRUNC, fcpl_id, fapl_id);
        if (file_id < 0) {
            H5Pclose(fapl_id);
            ESIO_ERROR("Unable to create file", ESIO_EFAILED);
        }
    } else {
        // Avoid overwriting existing file
        file_id = H5Fcreate(file, H5F_ACC_EXCL, fcpl_id, fapl_id);
        if (file_id < 0) {
            H5Pclose(fapl_id);
            ESIO_ERROR("File already exists", ESIO_EFAILED);
        }
    }

    // File creation successful: update state
    s->file_id = file_id;

    // Clean up temporary resources
    H5Pclose(fapl_id);

    return ESIO_SUCCESS;
}

int
esio_file_open(esio_state s, const char *file, int readwrite)
{
    // Sanity check incoming arguments
    if (s == NULL) {
        ESIO_ERROR("s == NULL", ESIO_EINVAL);
    }
    if (s->file_id != -1) {
        ESIO_ERROR("Cannot open new file because previous file not closed",
                   ESIO_EINVAL);
    }
    if (file == NULL) {
        ESIO_ERROR("file == NULL", ESIO_EINVAL);
    }

    // Initialize file access list property identifier
    const hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    if (fapl_id == -1) {
        ESIO_ERROR("Unable to create fapl_id", ESIO_ESANITY);
    }
    // Set collective details
    if (H5Pset_fapl_mpio(fapl_id, s->comm, s->info)) {
        H5Pclose(fapl_id);
        ESIO_ERROR("Unable to store MPI details in fapl_id", ESIO_ESANITY);
    }

    // Initialize access flags
    const unsigned flags = readwrite ? H5F_ACC_RDWR : H5F_ACC_RDONLY;

    // Collectively open the file
    const hid_t file_id = H5Fopen(file, flags, fapl_id);
    if (file_id < 0) {
        H5Pclose(fapl_id);
        ESIO_ERROR("Unable to open existing file", ESIO_EFAILED);
    }

    // File creation successful: update state
    s->file_id = file_id;

    // Clean up temporary resources
    H5Pclose(fapl_id);

    return ESIO_SUCCESS;
}

int esio_file_close(esio_state s)
{
    // Sanity check incoming arguments
    if (s == NULL) {
        ESIO_ERROR("s == NULL", ESIO_EINVAL);
    }
    if (s->file_id == -1) {
        ESIO_ERROR("No file currently open", ESIO_EINVAL);
    }

    // Close any open file
    if (H5Fclose(s->file_id) < 0) {
        ESIO_ERROR("Unable to close file", ESIO_EFAILED);
    }

    // File closing successful: update state
    s->file_id = -1;

    return ESIO_SUCCESS;
}

static
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

static
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

static
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
    return H5LTset_attribute_int(loc_id, name, "esio_metadata",
                                 metadata, ESIO_FIELD_METADATA_SIZE);
}

static
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
            loc_id, name, "esio_metadata", metadata);

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
        if (metadata[3] < 0 || metadata[3] >= esio_nlayout) {
            ESIO_ERROR_VAL("ESIO metadata contains unknown layout_tag",
                           ESIO_ESANITY, err);
        }
    }

    // Return the H5LTget_attribute_int error code
    return err;
}

static
hid_t esio_field_create(const esio_state s,
                        int cglobal, int bglobal, int aglobal,
                        const char* name, hid_t type_id)
{
    // Sanity check that the state's layout_tag matches our internal table
    if (esio_layout[s->layout_tag].tag != s->layout_tag) {
        ESIO_ERROR("SEVERE: Consistency error in esio_layout", ESIO_ESANITY);
    }

    // Create the filespace using current layout within state
    const hid_t filespace
        = (esio_layout[s->layout_tag].filespace_creator)(cglobal,
                                                         bglobal,
                                                         aglobal);
    if (filespace < 0) {
        ESIO_ERROR("Unable to create filespace", ESIO_ESANITY);
    }

    // Create the dataspace
    const hid_t dset_id
        = H5Dcreate1(s->file_id, name, type_id, filespace, H5P_DEFAULT);
    if (dset_id < 0) {
        ESIO_ERROR("Unable to create dataspace", ESIO_ESANITY);
    }

    // Stash field's metadata
    const herr_t status = esio_field_metadata_write(s->file_id, name,
                                                    s->layout_tag,
                                                    cglobal, bglobal, aglobal,
                                                    type_id);
    if (status < 0) {
        ESIO_ERROR("Unable to save field's ESIO metadata", ESIO_EFAILED);
    }

    // Clean up temporary resources
    H5Sclose(filespace);

    return dset_id;
}

static
int esio_field_close(hid_t dataset_id)
{
    if (H5Dclose(dataset_id) < 0) {
        ESIO_ERROR("Error closing field", ESIO_EFAILED);
    }

    return ESIO_SUCCESS;
}

int esio_field_size(const esio_state s,
                    const char* name,
                    int *cglobal, int *bglobal, int *aglobal)
{
    int ncomponents;
    const int status
        = esio_vfield_size(s, name, cglobal, bglobal, aglobal, &ncomponents);
    if (ncomponents != 1) {
        ESIO_ERROR("Named location must be treated as a vfield", ESIO_EINVAL);
    }
    return status;
}

int esio_vfield_size(const esio_state s,
                     const char* name,
                     int *cglobal, int *bglobal, int *aglobal,
                     int *ncomponents)
{
    // Sanity check incoming arguments
    if (s == NULL)        ESIO_ERROR("s == NULL",              ESIO_EINVAL);
    if (s->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EINVAL);

    const herr_t status = esio_field_metadata_read(
            s->file_id, name, NULL, cglobal, bglobal, aglobal, ncomponents);
    if (status < 0) {
        ESIO_ERROR("Unable to open field's ESIO metadata", ESIO_EFAILED);
    }

    return ESIO_SUCCESS;
}

static
int esio_field_write_internal(const esio_state s,
                              const char* name,
                              const void *field,
                              int cglobal, int cstart, int clocal, int cstride,
                              int bglobal, int bstart, int blocal, int bstride,
                              int aglobal, int astart, int alocal, int astride,
                              hid_t type_id)
{
    // Sanity check incoming arguments
    // Strides must be nonnegative because hsize_t is unsigned
    if (s == NULL)        ESIO_ERROR("s == NULL",              ESIO_EINVAL);
    if (s->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EINVAL);
    if (field == NULL)    ESIO_ERROR("field == NULL",          ESIO_EINVAL);
    if (cglobal  < 0)     ESIO_ERROR("cglobal < 0",            ESIO_EINVAL);
    if (cstart < 0)       ESIO_ERROR("cstart < 0",             ESIO_EINVAL);
    if (clocal < 1)       ESIO_ERROR("clocal < 1",             ESIO_EINVAL);
    if (cstride < 0)      ESIO_ERROR("cstride < 0",            ESIO_EINVAL);
    if (bglobal  < 0)     ESIO_ERROR("bglobal < 0",            ESIO_EINVAL);
    if (bstart < 0)       ESIO_ERROR("bstart < 0",             ESIO_EINVAL);
    if (blocal < 1)       ESIO_ERROR("blocal < 1",             ESIO_EINVAL);
    if (bstride < 0)      ESIO_ERROR("bstride < 0",            ESIO_EINVAL);
    if (aglobal  < 0)     ESIO_ERROR("aglobal < 0",            ESIO_EINVAL);
    if (astart < 0)       ESIO_ERROR("astart < 0",             ESIO_EINVAL);
    if (alocal < 1)       ESIO_ERROR("alocal < 1",             ESIO_EINVAL);
    if (astride < 0)      ESIO_ERROR("astride < 0",            ESIO_EINVAL);
    if (type_id < 0)      ESIO_ERROR("type_id < 0",            ESIO_EINVAL);

    // Provide contiguous defaults whenever the user supplied zero strides.
    // Strides are given in units of type_id; hence astride = 1 is contiguous.
    if (astride == 0) astride = 1;
    if (bstride == 0) bstride = astride * alocal;
    if (cstride == 0) cstride = bstride * blocal;

    // Attempt to read metadata for the field (which may or may not exist)
    int layout_tag;
    int field_cglobal, field_bglobal, field_aglobal;
    int field_ncomponents;
    const herr_t mstat = esio_field_metadata_read(s->file_id, name,
                                                  &layout_tag,
                                                  &field_cglobal,
                                                  &field_bglobal,
                                                  &field_aglobal,
                                                  &field_ncomponents);

    if (mstat < 0) {
        // Field did not exist

        // Create dataset and write it with the active field layout
        const hid_t dset_id
            = esio_field_create(s, cglobal, bglobal, aglobal, name, type_id);
        const int wstat = (esio_layout[s->layout_tag].field_writer)(
                dset_id, field,
                cglobal, cstart, clocal, cstride,
                bglobal, bstart, blocal, bstride,
                aglobal, astart, alocal, astride,
                type_id);
        if (wstat != ESIO_SUCCESS) {
            esio_field_close(dset_id);
            ESIO_ERROR_VAL("Error writing new field", ESIO_EFAILED, wstat);
        }
        esio_field_close(dset_id);

    } else {
        // Field already existed

        // Ensure caller gave correct size information
        if (cglobal != field_cglobal) {
            ESIO_ERROR("request cglobal mismatch with existing field",
                       ESIO_EINVAL);
        }
        if (bglobal != field_bglobal) {
            ESIO_ERROR("request bglobal mismatch with existing field",
                       ESIO_EINVAL);
        }
        if (aglobal != field_aglobal) {
            ESIO_ERROR("request aglobal mismatch with existing field",
                       ESIO_EINVAL);
        }

        // Ensure caller gave type with correct component count
        if (esio_type_ncomponents(type_id) != field_ncomponents) {
            ESIO_ERROR("request ncomponents mismatch with existing field",
                       ESIO_EINVAL);
        }

        // Open the existing field's dataset
        const hid_t dset_id = H5Dopen1(s->file_id, name);
        if (dset_id < 0) {
            ESIO_ERROR("Unable to open dataset", ESIO_EFAILED);
        }

        // Check if supplied type can be converted to the field's type
        const hid_t field_type_id = H5Dget_type(dset_id);
        H5T_cdata_t *pcdata;
        const H5T_conv_t converter = H5Tfind(type_id, field_type_id, &pcdata);
        if (converter == NULL) {
            H5Tclose(field_type_id);
            ESIO_ERROR("request type not convertible to existing field type",
                       ESIO_EINVAL);
        }
        H5Tclose(field_type_id);

        // Overwrite existing data using layout routines per metadata
        const int wstat = (esio_layout[layout_tag].field_writer)(
                dset_id, field,
                cglobal, cstart, clocal, cstride,
                bglobal, bstart, blocal, bstride,
                aglobal, astart, alocal, astride,
                type_id);
        if (wstat != ESIO_SUCCESS) {
            esio_field_close(dset_id);
            ESIO_ERROR_VAL("Error writing new field", ESIO_EFAILED, wstat);
        }
        esio_field_close(dset_id);

    }

    return ESIO_SUCCESS;
}

static
int esio_field_read_internal(const esio_state s,
                             const char* name,
                             void *field,
                             int cglobal, int cstart, int clocal, int cstride,
                             int bglobal, int bstart, int blocal, int bstride,
                             int aglobal, int astart, int alocal, int astride,
                             hid_t type_id)
{
    // Sanity check incoming arguments
    // Strides must be nonnegative because hsize_t is unsigned
    if (s == NULL)        ESIO_ERROR("s == NULL",              ESIO_EINVAL);
    if (s->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EINVAL);
    if (field == NULL)    ESIO_ERROR("field == NULL",          ESIO_EINVAL);
    if (cglobal  < 0)     ESIO_ERROR("cglobal < 0",            ESIO_EINVAL);
    if (cstart < 0)       ESIO_ERROR("cstart < 0",             ESIO_EINVAL);
    if (clocal < 1)       ESIO_ERROR("clocal < 1",             ESIO_EINVAL);
    if (cstride < 0)      ESIO_ERROR("cstride < 0",            ESIO_EINVAL);
    if (bglobal  < 0)     ESIO_ERROR("bglobal < 0",            ESIO_EINVAL);
    if (bstart < 0)       ESIO_ERROR("bstart < 0",             ESIO_EINVAL);
    if (blocal < 1)       ESIO_ERROR("blocal < 1",             ESIO_EINVAL);
    if (bstride < 0)      ESIO_ERROR("bstride < 0",            ESIO_EINVAL);
    if (aglobal  < 0)     ESIO_ERROR("aglobal < 0",            ESIO_EINVAL);
    if (astart < 0)       ESIO_ERROR("astart < 0",             ESIO_EINVAL);
    if (alocal < 1)       ESIO_ERROR("alocal < 1",             ESIO_EINVAL);
    if (astride < 0)      ESIO_ERROR("astride < 0",            ESIO_EINVAL);
    if (type_id < 0)      ESIO_ERROR("type_id < 0",            ESIO_EINVAL);

    // Provide contiguous defaults whenever the user supplied zero strides.
    // Strides are given in units of type_id; hence astride = 1 is contiguous.
    if (astride == 0) astride = 1;
    if (bstride == 0) bstride = astride * alocal;
    if (cstride == 0) cstride = bstride * blocal;

    // Read metadata for the field
    int layout_tag;
    int field_cglobal, field_bglobal, field_aglobal;
    int field_ncomponents;
    const herr_t status = esio_field_metadata_read(s->file_id, name,
                                                   &layout_tag,
                                                   &field_cglobal,
                                                   &field_bglobal,
                                                   &field_aglobal,
                                                   &field_ncomponents);
    if (status < 0) {
        ESIO_ERROR("Unable to read field's ESIO metadata", ESIO_EFAILED);
    }

    // Ensure caller gave correct size information
    if (cglobal != field_cglobal) {
        ESIO_ERROR("field read request has incorrect cglobal", ESIO_EINVAL);
    }
    if (bglobal != field_bglobal) {
        ESIO_ERROR("field read request has incorrect bglobal", ESIO_EINVAL);
    }
    if (aglobal != field_aglobal) {
        ESIO_ERROR("field read request has incorrect aglobal", ESIO_EINVAL);
    }

    // Ensure caller gave type with correct component count
    if (esio_type_ncomponents(type_id) != field_ncomponents) {
        ESIO_ERROR("request ncomponents mismatch with existing field",
                    ESIO_EINVAL);
    }

    // Open existing dataset
    const hid_t dapl_id = H5P_DEFAULT;
    const hid_t dset_id = H5Dopen2(s->file_id, name, dapl_id);
    if (dset_id < 0) {
        ESIO_ERROR("Unable to open dataset", ESIO_EFAILED);
    }

    // Check if supplied type can be converted to the field's type
    const hid_t field_type_id = H5Dget_type(dset_id);
    H5T_cdata_t *pcdata;
    const H5T_conv_t converter = H5Tfind(type_id, field_type_id, &pcdata);
    if (converter == NULL) {
        H5Tclose(field_type_id);
        ESIO_ERROR("request type not convertible to existing field type",
                    ESIO_EINVAL);
    }
    H5Tclose(field_type_id);

    // Read the field based on the metadata's layout_tag
    // Note that this means we can read any layout ESIO understands
    // Note that reading does not change the chosen field write layout_tag
    (esio_layout[layout_tag].field_reader)(dset_id, field,
                                           cglobal, cstart, clocal, cstride,
                                           bglobal, bstart, blocal, bstride,
                                           aglobal, astart, alocal, astride,
                                           type_id);

    // Close dataset
    esio_field_close(dset_id);

    return ESIO_SUCCESS;
}

#define GEN_FIELD_OP(OP,QUAL,TYPE,H5TYPE)                                   \
int esio_field_ ## OP ## _ ## TYPE (                                        \
        const esio_state s,                                                 \
        const char* name,                                                   \
        QUAL TYPE *field,                                                   \
        int cglobal, int cstart, int clocal, int cstride,                   \
        int bglobal, int bstart, int blocal, int bstride,                   \
        int aglobal, int astart, int alocal, int astride)                   \
{                                                                           \
    return esio_field_ ## OP ## _internal(s, name, field,                   \
                                          cglobal, cstart, clocal, cstride, \
                                          bglobal, bstart, blocal, bstride, \
                                          aglobal, astart, alocal, astride, \
                                          H5TYPE);                          \
}

GEN_FIELD_OP(write, const,       double, H5T_NATIVE_DOUBLE)
GEN_FIELD_OP(read,  /*mutable*/, double, H5T_NATIVE_DOUBLE)

GEN_FIELD_OP(write, const,       float, H5T_NATIVE_FLOAT)
GEN_FIELD_OP(read,  /*mutable*/, float, H5T_NATIVE_FLOAT)


#define GEN_VFIELD_OP(OP,QUAL,TYPE,H5TYPE)                                \
int esio_vfield_ ## OP ## _ ## TYPE(                                      \
        const esio_state s,                                               \
        const char* name,                                                 \
        QUAL TYPE *field,                                                 \
        int cglobal, int cstart, int clocal, int cstride,                 \
        int bglobal, int bstart, int blocal, int bstride,                 \
        int aglobal, int astart, int alocal, int astride,                 \
        int ncomponents)                                                  \
{                                                                         \
    const hid_t array_type_id = esio_type_arrayify(H5TYPE, ncomponents);  \
                                                                          \
    /* Input strides are in units of sizeof(TYPE). HDF5 requires slab */  \
    /* selection using sizeof(array_type_id) and cannot use "partial" */  \
    /* offsets in this selection process.  Check strides conform.     */  \
    if (cstride % ncomponents) {                                          \
        H5Tclose(array_type_id);                                          \
        ESIO_ERROR("cstride must be an integer multiple of ncomponents",  \
                   ESIO_EINVAL);                                          \
    }                                                                     \
    if (bstride % ncomponents) {                                          \
        H5Tclose(array_type_id);                                          \
        ESIO_ERROR("bstride must be an integer multiple of ncomponents",  \
                   ESIO_EINVAL);                                          \
    }                                                                     \
    if (astride % ncomponents) {                                          \
        H5Tclose(array_type_id);                                          \
        ESIO_ERROR("astride must be an integer multiple of ncomponents",  \
                   ESIO_EINVAL);                                          \
    }                                                                     \
    const int retval = esio_field_ ## OP ## _internal(                    \
            s, name, field,                                               \
            cglobal, cstart, clocal, (cstride / ncomponents),             \
            bglobal, bstart, blocal, (bstride / ncomponents),             \
            aglobal, astart, alocal, (astride / ncomponents),             \
            array_type_id);                                               \
    H5Tclose(array_type_id);                                              \
    return retval;                                                        \
}

GEN_VFIELD_OP(write, const,       double, H5T_NATIVE_DOUBLE)
GEN_VFIELD_OP(read,  /*mutable*/, double, H5T_NATIVE_DOUBLE)

GEN_VFIELD_OP(write, const,       float, H5T_NATIVE_FLOAT)
GEN_VFIELD_OP(read,  /*mutable*/, float, H5T_NATIVE_FLOAT)
