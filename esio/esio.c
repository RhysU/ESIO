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
#include "metadata.h"
#include "version.h"

//*********************************************************************
// INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL
//*********************************************************************

static
MPI_Comm esio_MPI_Comm_dup_with_name(MPI_Comm comm);

static
hid_t esio_field_create(const esio_state s,
                        int cglobal, int bglobal, int aglobal,
                        const char* name, hid_t type_id);

static
hid_t esio_plane_create(const esio_state s,
                        int bglobal, int aglobal,
                        const char* name, hid_t type_id);

static
hid_t esio_line_create(const esio_state s,
                       int aglobal,
                       const char* name, hid_t type_id);

static
int esio_field_close(hid_t dataset_id);

static
int esio_plane_close(hid_t dataset_id);

static
int esio_line_close(hid_t dataset_id);

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

static
int esio_plane_write_internal(const esio_state s,
                              const char* name,
                              const void *plane,
                              int bglobal, int bstart, int blocal, int bstride,
                              int aglobal, int astart, int alocal, int astride,
                              hid_t type_id);

static
int esio_plane_read_internal(const esio_state s,
                             const char* name,
                             void *plane,
                             int bglobal, int bstart, int blocal, int bstride,
                             int aglobal, int astart, int alocal, int astride,
                             hid_t type_id);

static
int esio_line_write_internal(const esio_state s,
                             const char* name,
                             const void *line,
                             int aglobal, int astart, int alocal, int astride,
                             hid_t type_id);

static
int esio_line_read_internal(const esio_state s,
                            const char* name,
                            void *line,
                            int aglobal, int astart, int alocal, int astride,
                            hid_t type_id);

static
int esio_point_write_internal(const esio_state s,
                              const char* name,
                              const void *point,
                              hid_t type_id);

static
int esio_point_read_internal(const esio_state s,
                             const char* name,
                             void *point,
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
    },
    {
        1,
        &esio_layout1_filespace_creator,
        &esio_layout1_field_writer,
        &esio_layout1_field_reader
    },
};
const int esio_nlayout = sizeof(esio_layout)/sizeof(esio_layout[0]);


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

int esio_layout_count()
{
    return esio_nlayout;
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
hid_t esio_plane_create(const esio_state s,
                        int bglobal, int aglobal,
                        const char* name, hid_t type_id)
{
    // Create the filespace
    const hsize_t dims[2] = { bglobal, aglobal };
    const hid_t filespace = H5Screate_simple(2, dims, NULL);
    if (filespace < 0) {
        ESIO_ERROR("Unable to create filespace", ESIO_ESANITY);
    }

    // Create the dataspace
    const hid_t dset_id
        = H5Dcreate1(s->file_id, name, type_id, filespace, H5P_DEFAULT);
    if (dset_id < 0) {
        ESIO_ERROR("Unable to create dataspace", ESIO_ESANITY);
    }

    // Stash plane's metadata
    const herr_t status = esio_plane_metadata_write(s->file_id, name,
                                                    bglobal, aglobal,
                                                    type_id);
    if (status < 0) {
        ESIO_ERROR("Unable to save plane's ESIO metadata", ESIO_EFAILED);
    }

    // Clean up temporary resources
    H5Sclose(filespace);

    return dset_id;
}

static
hid_t esio_line_create(const esio_state s,
                       int aglobal,
                       const char* name, hid_t type_id)
{
    // Create the filespace
    const hsize_t dims[1] = { aglobal };
    const hid_t filespace = H5Screate_simple(1, dims, NULL);
    if (filespace < 0) {
        ESIO_ERROR("Unable to create filespace", ESIO_ESANITY);
    }

    // Create the dataspace
    const hid_t dset_id
        = H5Dcreate1(s->file_id, name, type_id, filespace, H5P_DEFAULT);
    if (dset_id < 0) {
        ESIO_ERROR("Unable to create dataspace", ESIO_ESANITY);
    }

    // Stash line's metadata
    const herr_t status = esio_line_metadata_write(s->file_id, name,
                                                   aglobal,
                                                   type_id);
    if (status < 0) {
        ESIO_ERROR("Unable to save line's ESIO metadata", ESIO_EFAILED);
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

static
int esio_plane_close(hid_t dataset_id)
{
    if (H5Dclose(dataset_id) < 0) {
        ESIO_ERROR("Error closing plane", ESIO_EFAILED);
    }

    return ESIO_SUCCESS;
}

static
int esio_line_close(hid_t dataset_id)
{
    if (H5Dclose(dataset_id) < 0) {
        ESIO_ERROR("Error closing line", ESIO_EFAILED);
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

int esio_plane_size(const esio_state s,
                    const char* name,
                    int *bglobal, int *aglobal)
{
    int ncomponents;
    const int status
        = esio_vplane_size(s, name, bglobal, aglobal, &ncomponents);
    if (ncomponents != 1) {
        ESIO_ERROR("Named location must be treated as a vplane", ESIO_EINVAL);
    }
    return status;
}

int esio_vplane_size(const esio_state s,
                     const char* name,
                     int *bglobal, int *aglobal,
                     int *ncomponents)
{
    // Sanity check incoming arguments
    if (s == NULL)        ESIO_ERROR("s == NULL",              ESIO_EINVAL);
    if (s->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EINVAL);

    const herr_t status = esio_plane_metadata_read(
            s->file_id, name, bglobal, aglobal, ncomponents);
    if (status < 0) {
        ESIO_ERROR("Unable to open plane's ESIO metadata", ESIO_EFAILED);
    }

    return ESIO_SUCCESS;
}

int esio_line_size(const esio_state s,
                    const char* name,
                    int *aglobal)
{
    int ncomponents;
    const int status = esio_vline_size(s, name, aglobal, &ncomponents);
    if (ncomponents != 1) {
        ESIO_ERROR("Named location must be treated as a vline", ESIO_EINVAL);
    }
    return status;
}

int esio_vline_size(const esio_state s,
                    const char* name,
                    int *aglobal,
                    int *ncomponents)
{
    // Sanity check incoming arguments
    if (s == NULL)        ESIO_ERROR("s == NULL",              ESIO_EINVAL);
    if (s->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EINVAL);

    const herr_t status = esio_line_metadata_read(
            s->file_id, name, aglobal, ncomponents);
    if (status < 0) {
        ESIO_ERROR("Unable to open line's ESIO metadata", ESIO_EFAILED);
    }

    return ESIO_SUCCESS;
}

// *******************************************************************
// FIELD READ WRITE FIELD READ WRITE FIELD READ WRITE FIELD READ WRITE
// *******************************************************************

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
            H5Dclose(dset_id);
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
            ESIO_ERROR_VAL("Error writing overwriting field", ESIO_EFAILED, wstat);
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
        H5Dclose(dset_id);
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

// *******************************************************************
// PLANE READ WRITE PLANE READ WRITE PLANE READ WRITE PLANE READ WRITE
// *******************************************************************

static
int esio_plane_write_internal(const esio_state s,
                              const char* name,
                              const void *plane,
                              int bglobal, int bstart, int blocal, int bstride,
                              int aglobal, int astart, int alocal, int astride,
                              hid_t type_id)
{
    // Sanity check incoming arguments
    // Strides must be nonnegative because hsize_t is unsigned
    if (s == NULL)        ESIO_ERROR("s == NULL",              ESIO_EINVAL);
    if (s->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EINVAL);
    if (plane == NULL)    ESIO_ERROR("plane == NULL",          ESIO_EINVAL);
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

    // Attempt to read metadata for the plane (which may or may not exist)
    int plane_bglobal, plane_aglobal;
    int plane_ncomponents;
    const herr_t mstat = esio_plane_metadata_read(s->file_id, name,
                                                  &plane_bglobal,
                                                  &plane_aglobal,
                                                  &plane_ncomponents);

    hid_t dset_id;
    if (mstat < 0) {
        // Plane did not exist so create it
        dset_id = esio_plane_create(s, bglobal, aglobal, name, type_id);
        if (dset_id < 0) {
            ESIO_ERROR("Error creating new plane", ESIO_EFAILED);
        }
    } else {
        // Plane already existed

        // Ensure caller gave correct size information
        if (bglobal != plane_bglobal) {
            ESIO_ERROR("request bglobal mismatch with existing plane",
                       ESIO_EINVAL);
        }
        if (aglobal != plane_aglobal) {
            ESIO_ERROR("request aglobal mismatch with existing plane",
                       ESIO_EINVAL);
        }

        // Ensure caller gave type with correct component count
        if (esio_type_ncomponents(type_id) != plane_ncomponents) {
            ESIO_ERROR("request ncomponents mismatch with existing plane",
                       ESIO_EINVAL);
        }

        // Open the existing plane's dataset
        dset_id = H5Dopen1(s->file_id, name);
        if (dset_id < 0) {
            ESIO_ERROR("Unable to open dataset", ESIO_EFAILED);
        }

        // Check if supplied type can be converted to the plane's type
        const hid_t plane_type_id = H5Dget_type(dset_id);
        H5T_cdata_t *pcdata;
        const H5T_conv_t converter = H5Tfind(type_id, plane_type_id, &pcdata);
        if (converter == NULL) {
            H5Tclose(plane_type_id);
            H5Dclose(dset_id);
            ESIO_ERROR("request type not convertible to existing plane type",
                       ESIO_EINVAL);
        }
        H5Tclose(plane_type_id);
    }

    // Write field
    const int wstat = esio_plane_writer(dset_id, plane,
                                        bglobal, bstart, blocal, bstride,
                                        aglobal, astart, alocal, astride,
                                        type_id);
    if (wstat != ESIO_SUCCESS) {
        esio_plane_close(dset_id);
        ESIO_ERROR_VAL("Error writing plane", ESIO_EFAILED, wstat);
    }
    esio_plane_close(dset_id);

    return ESIO_SUCCESS;
}

static
int esio_plane_read_internal(const esio_state s,
                             const char* name,
                             void *plane,
                             int bglobal, int bstart, int blocal, int bstride,
                             int aglobal, int astart, int alocal, int astride,
                             hid_t type_id)
{
    // Sanity check incoming arguments
    // Strides must be nonnegative because hsize_t is unsigned
    if (s == NULL)        ESIO_ERROR("s == NULL",              ESIO_EINVAL);
    if (s->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EINVAL);
    if (plane == NULL)    ESIO_ERROR("plane == NULL",          ESIO_EINVAL);
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

    // Read metadata for the plane
    int plane_bglobal, plane_aglobal;
    int plane_ncomponents;
    const herr_t status = esio_plane_metadata_read(s->file_id, name,
                                                   &plane_bglobal,
                                                   &plane_aglobal,
                                                   &plane_ncomponents);
    if (status < 0) {
        ESIO_ERROR("Unable to read plane's ESIO metadata", ESIO_EFAILED);
    }

    // Ensure caller gave correct size information
    if (bglobal != plane_bglobal) {
        ESIO_ERROR("plane read request has incorrect bglobal", ESIO_EINVAL);
    }
    if (aglobal != plane_aglobal) {
        ESIO_ERROR("plane read request has incorrect aglobal", ESIO_EINVAL);
    }

    // Ensure caller gave type with correct component count
    if (esio_type_ncomponents(type_id) != plane_ncomponents) {
        ESIO_ERROR("request ncomponents mismatch with existing plane",
                    ESIO_EINVAL);
    }

    // Open existing dataset
    const hid_t dapl_id = H5P_DEFAULT;
    const hid_t dset_id = H5Dopen2(s->file_id, name, dapl_id);
    if (dset_id < 0) {
        ESIO_ERROR("Unable to open dataset", ESIO_EFAILED);
    }

    // Check if supplied type can be converted to the plane's type
    const hid_t plane_type_id = H5Dget_type(dset_id);
    H5T_cdata_t *pcdata;
    const H5T_conv_t converter = H5Tfind(type_id, plane_type_id, &pcdata);
    if (converter == NULL) {
        H5Tclose(plane_type_id);
        H5Dclose(dset_id);
        ESIO_ERROR("request type not convertible to existing plane type",
                    ESIO_EINVAL);
    }
    H5Tclose(plane_type_id);

    // Read plane
    const int rstat = esio_plane_reader(dset_id, plane,
                                        bglobal, bstart, blocal, bstride,
                                        aglobal, astart, alocal, astride,
                                        type_id);
    if (rstat != ESIO_SUCCESS) {
        esio_plane_close(dset_id);
        ESIO_ERROR_VAL("Error reading plane", ESIO_EFAILED, rstat);
    }
    esio_plane_close(dset_id);

    return ESIO_SUCCESS;
}

#define GEN_PLANE_OP(OP,QUAL,TYPE,H5TYPE)                                   \
int esio_plane_ ## OP ## _ ## TYPE (                                        \
        const esio_state s,                                                 \
        const char* name,                                                   \
        QUAL TYPE *plane,                                                   \
        int bglobal, int bstart, int blocal, int bstride,                   \
        int aglobal, int astart, int alocal, int astride)                   \
{                                                                           \
    return esio_plane_ ## OP ## _internal(s, name, plane,                   \
                                          bglobal, bstart, blocal, bstride, \
                                          aglobal, astart, alocal, astride, \
                                          H5TYPE);                          \
}

GEN_PLANE_OP(write, const,       double, H5T_NATIVE_DOUBLE)
GEN_PLANE_OP(read,  /*mutable*/, double, H5T_NATIVE_DOUBLE)

GEN_PLANE_OP(write, const,       float, H5T_NATIVE_FLOAT)
GEN_PLANE_OP(read,  /*mutable*/, float, H5T_NATIVE_FLOAT)


#define GEN_VPLANE_OP(OP,QUAL,TYPE,H5TYPE)                                \
int esio_vplane_ ## OP ## _ ## TYPE(                                      \
        const esio_state s,                                               \
        const char* name,                                                 \
        QUAL TYPE *plane,                                                 \
        int bglobal, int bstart, int blocal, int bstride,                 \
        int aglobal, int astart, int alocal, int astride,                 \
        int ncomponents)                                                  \
{                                                                         \
    const hid_t array_type_id = esio_type_arrayify(H5TYPE, ncomponents);  \
                                                                          \
    /* Input strides are in units of sizeof(TYPE). HDF5 requires slab */  \
    /* selection using sizeof(array_type_id) and cannot use "partial" */  \
    /* offsets in this selection process.  Check strides conform.     */  \
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
    const int retval = esio_plane_ ## OP ## _internal(                    \
            s, name, plane,                                               \
            bglobal, bstart, blocal, (bstride / ncomponents),             \
            aglobal, astart, alocal, (astride / ncomponents),             \
            array_type_id);                                               \
    H5Tclose(array_type_id);                                              \
    return retval;                                                        \
}

GEN_VPLANE_OP(write, const,       double, H5T_NATIVE_DOUBLE)
GEN_VPLANE_OP(read,  /*mutable*/, double, H5T_NATIVE_DOUBLE)

GEN_VPLANE_OP(write, const,       float, H5T_NATIVE_FLOAT)
GEN_VPLANE_OP(read,  /*mutable*/, float, H5T_NATIVE_FLOAT)

// *******************************************************************
// LINE READ WRITE LINE READ WRITE LINE READ WRITE LINE READ WRITE
// *******************************************************************

static
int esio_line_write_internal(const esio_state s,
                             const char* name,
                             const void *line,
                             int aglobal, int astart, int alocal, int astride,
                             hid_t type_id)
{
    // Sanity check incoming arguments
    // Strides must be nonnegative because hsize_t is unsigned
    if (s == NULL)        ESIO_ERROR("s == NULL",              ESIO_EINVAL);
    if (s->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EINVAL);
    if (line == NULL)     ESIO_ERROR("line == NULL",           ESIO_EINVAL);
    if (aglobal  < 0)     ESIO_ERROR("aglobal < 0",            ESIO_EINVAL);
    if (astart < 0)       ESIO_ERROR("astart < 0",             ESIO_EINVAL);
    if (alocal < 1)       ESIO_ERROR("alocal < 1",             ESIO_EINVAL);
    if (astride < 0)      ESIO_ERROR("astride < 0",            ESIO_EINVAL);
    if (type_id < 0)      ESIO_ERROR("type_id < 0",            ESIO_EINVAL);

    // Provide contiguous defaults whenever the user supplied zero strides.
    // Strides are given in units of type_id; hence astride = 1 is contiguous.
    if (astride == 0) astride = 1;

    // Attempt to read metadata for the line (which may or may not exist)
    int line_aglobal, line_ncomponents;
    const herr_t mstat = esio_line_metadata_read(s->file_id, name,
                                                 &line_aglobal,
                                                 &line_ncomponents);

    hid_t dset_id;
    if (mstat < 0) {
        // Line did not exist so create it
        dset_id = esio_line_create(s, aglobal, name, type_id);
        if (dset_id < 0) {
            ESIO_ERROR("Error creating new line", ESIO_EFAILED);
        }
    } else {
        // Line already existed

        // Ensure caller gave correct size information
        if (aglobal != line_aglobal) {
            ESIO_ERROR("request aglobal mismatch with existing line",
                       ESIO_EINVAL);
        }

        // Ensure caller gave type with correct component count
        if (esio_type_ncomponents(type_id) != line_ncomponents) {
            ESIO_ERROR("request ncomponents mismatch with existing line",
                       ESIO_EINVAL);
        }

        // Open the existing line's dataset
        dset_id = H5Dopen1(s->file_id, name);
        if (dset_id < 0) {
            ESIO_ERROR("Unable to open dataset", ESIO_EFAILED);
        }

        // Check if supplied type can be converted to the line's type
        const hid_t line_type_id = H5Dget_type(dset_id);
        H5T_cdata_t *pcdata;
        const H5T_conv_t converter = H5Tfind(type_id, line_type_id, &pcdata);
        if (converter == NULL) {
            H5Tclose(line_type_id);
            H5Dclose(dset_id);
            ESIO_ERROR("request type not convertible to existing line type",
                       ESIO_EINVAL);
        }
        H5Tclose(line_type_id);
    }

    // Write field
    const int wstat = esio_line_writer(dset_id, line,
                                       aglobal, astart, alocal, astride,
                                       type_id);
    if (wstat != ESIO_SUCCESS) {
        esio_line_close(dset_id);
        ESIO_ERROR_VAL("Error writing line", ESIO_EFAILED, wstat);
    }
    esio_line_close(dset_id);

    return ESIO_SUCCESS;
}

static
int esio_line_read_internal(const esio_state s,
                            const char* name,
                            void *line,
                            int aglobal, int astart, int alocal, int astride,
                            hid_t type_id)
{
    // Sanity check incoming arguments
    // Strides must be nonnegative because hsize_t is unsigned
    if (s == NULL)        ESIO_ERROR("s == NULL",              ESIO_EINVAL);
    if (s->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EINVAL);
    if (line == NULL)     ESIO_ERROR("line == NULL",           ESIO_EINVAL);
    if (aglobal  < 0)     ESIO_ERROR("aglobal < 0",            ESIO_EINVAL);
    if (astart < 0)       ESIO_ERROR("astart < 0",             ESIO_EINVAL);
    if (alocal < 1)       ESIO_ERROR("alocal < 1",             ESIO_EINVAL);
    if (astride < 0)      ESIO_ERROR("astride < 0",            ESIO_EINVAL);
    if (type_id < 0)      ESIO_ERROR("type_id < 0",            ESIO_EINVAL);

    // Provide contiguous defaults whenever the user supplied zero strides.
    // Strides are given in units of type_id; hence astride = 1 is contiguous.
    if (astride == 0) astride = 1;

    // Read metadata for the line
    int line_aglobal, line_ncomponents;
    const herr_t status = esio_line_metadata_read(s->file_id, name,
                                                  &line_aglobal,
                                                  &line_ncomponents);
    if (status < 0) {
        ESIO_ERROR("Unable to read line's ESIO metadata", ESIO_EFAILED);
    }

    // Ensure caller gave correct size information
    if (aglobal != line_aglobal) {
        ESIO_ERROR("line read request has incorrect aglobal", ESIO_EINVAL);
    }

    // Ensure caller gave type with correct component count
    if (esio_type_ncomponents(type_id) != line_ncomponents) {
        ESIO_ERROR("request ncomponents mismatch with existing line",
                    ESIO_EINVAL);
    }

    // Open existing dataset
    const hid_t dapl_id = H5P_DEFAULT;
    const hid_t dset_id = H5Dopen2(s->file_id, name, dapl_id);
    if (dset_id < 0) {
        ESIO_ERROR("Unable to open dataset", ESIO_EFAILED);
    }

    // Check if supplied type can be converted to the line's type
    const hid_t line_type_id = H5Dget_type(dset_id);
    H5T_cdata_t *pcdata;
    const H5T_conv_t converter = H5Tfind(type_id, line_type_id, &pcdata);
    if (converter == NULL) {
        H5Tclose(line_type_id);
        H5Dclose(dset_id);
        ESIO_ERROR("request type not convertible to existing line type",
                    ESIO_EINVAL);
    }
    H5Tclose(line_type_id);

    // Read line
    const int rstat = esio_line_reader(dset_id, line,
                                       aglobal, astart, alocal, astride,
                                       type_id);
    if (rstat != ESIO_SUCCESS) {
        esio_line_close(dset_id);
        ESIO_ERROR_VAL("Error reading line", ESIO_EFAILED, rstat);
    }
    esio_line_close(dset_id);

    return ESIO_SUCCESS;
}

#define GEN_LINE_OP(OP,QUAL,TYPE,H5TYPE)                                    \
int esio_line_ ## OP ## _ ## TYPE (                                         \
        const esio_state s,                                                 \
        const char* name,                                                   \
        QUAL TYPE *line,                                                    \
        int aglobal, int astart, int alocal, int astride)                   \
{                                                                           \
    return esio_line_ ## OP ## _internal(s, name, line,                     \
                                          aglobal, astart, alocal, astride, \
                                          H5TYPE);                          \
}

GEN_LINE_OP(write, const,       double, H5T_NATIVE_DOUBLE)
GEN_LINE_OP(read,  /*mutable*/, double, H5T_NATIVE_DOUBLE)

GEN_LINE_OP(write, const,       float, H5T_NATIVE_FLOAT)
GEN_LINE_OP(read,  /*mutable*/, float, H5T_NATIVE_FLOAT)


#define GEN_VLINE_OP(OP,QUAL,TYPE,H5TYPE)                                 \
int esio_vline_ ## OP ## _ ## TYPE(                                       \
        const esio_state s,                                               \
        const char* name,                                                 \
        QUAL TYPE *line,                                                  \
        int aglobal, int astart, int alocal, int astride,                 \
        int ncomponents)                                                  \
{                                                                         \
    const hid_t array_type_id = esio_type_arrayify(H5TYPE, ncomponents);  \
                                                                          \
    /* Input strides are in units of sizeof(TYPE). HDF5 requires slab */  \
    /* selection using sizeof(array_type_id) and cannot use "partial" */  \
    /* offsets in this selection process.  Check strides conform.     */  \
    if (astride % ncomponents) {                                          \
        H5Tclose(array_type_id);                                          \
        ESIO_ERROR("astride must be an integer multiple of ncomponents",  \
                   ESIO_EINVAL);                                          \
    }                                                                     \
    const int retval = esio_line_ ## OP ## _internal(                     \
            s, name, line,                                                \
            aglobal, astart, alocal, (astride / ncomponents),             \
            array_type_id);                                               \
    H5Tclose(array_type_id);                                              \
    return retval;                                                        \
}

GEN_VLINE_OP(write, const,       double, H5T_NATIVE_DOUBLE)
GEN_VLINE_OP(read,  /*mutable*/, double, H5T_NATIVE_DOUBLE)

GEN_VLINE_OP(write, const,       float, H5T_NATIVE_FLOAT)
GEN_VLINE_OP(read,  /*mutable*/, float, H5T_NATIVE_FLOAT)
