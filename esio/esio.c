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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <mpi.h>
#include "esio.h"
#include "error.h"

//*********************************************************************
// INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL
//*********************************************************************

static
MPI_Comm esio_MPI_Comm_dup_with_name(MPI_Comm comm);

static
hid_t layout1_filespace_creator(int na, int nb, int nc);

static
int layout1_field_writer(hid_t dset_id, const void *field,
                         int na, int ast, int asz,
                         int nb, int bst, int bsz,
                         int nc, int cst, int csz,
                         hid_t type_id, size_t type_size);

static
hid_t esio_field_create(esio_state s,
                        int na, int nb, int nc,
                        const char* name, hid_t type_id);

static
int esio_field_close(hid_t dataset_id);

static
int esio_field_write_internal(esio_state s,
                              const char* name,
                              const void *field,
                              int na, int ast, int asz,
                              int nb, int bst, int bsz,
                              int nc, int cst, int csz,
                              hid_t type_id,
                              size_t type_size);

//*********************************************************************
// INTERNAL TYPES INTERNAL TYPES INTERNAL TYPES INTERNAL TYPES INTERNAL
//*********************************************************************

typedef struct esio_layout {
    int tag;
    const char * name;
    hid_t (*filespace_creator)(int, int, int);
    int   (*field_writer)     (hid_t, const void *,
                               int, int, int,
                               int, int, int,
                               int, int, int,
                               hid_t, size_t);
} esio_layout;

static const esio_layout layout_1 = {
    1, "2D(a_bc)", &layout1_filespace_creator, &layout1_field_writer
};

struct esio_state_s {
    MPI_Comm    comm;      //< Communicator used for collective operations
    int         comm_rank; //< Process rank within in MPI communicator
    int         comm_size; //< Number of ranks within MPI communicator
    MPI_Info    info;      //< Info object used for collective operations
    hid_t       file_id;   //< Active HDF file identifier
    esio_layout layout;    //< Active field layout within HDF5 file
};


//***************************************************************************
// IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION
//***************************************************************************

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
    s->comm              = esio_MPI_Comm_dup_with_name(comm);
    s->comm_rank         = comm_rank;
    s->comm_size         = comm_size;
    s->info              = info;
    s->file_id           = -1;
    s->layout            = layout_1;

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
            ESIO_ERROR("Unable to create file", ESIO_ESANITY);
        }
    } else {
        // Avoid overwriting existing file
        file_id = H5Fcreate(file, H5F_ACC_EXCL, fcpl_id, fapl_id);
        if (file_id < 0) {
            H5Pclose(fapl_id);
            ESIO_ERROR("File already exists", ESIO_ESANITY);
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
        ESIO_ERROR("Unable to open existing file", ESIO_ESANITY);
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
hid_t esio_field_create(esio_state s,
                        int na, int nb, int nc,
                        const char* name, hid_t type_id)
{
    // Create the filespace using current layout within state state
    const hid_t filespace = (s->layout.filespace_creator)(na, nb, nc);
    if (filespace < 0) {
        ESIO_ERROR("Unable to create filespace", ESIO_ESANITY);
    }

    // Create the dataspace
    const hid_t dset_id
        = H5Dcreate1(s->file_id, name, type_id, filespace, H5P_DEFAULT);
    if (dset_id < 0) {
        ESIO_ERROR("Unable to create dataspace", ESIO_ESANITY);
    }

    // Stash the layout tag as an attribute
    if (H5LTset_attribute_int(s->file_id, name, "esio_layout_tag",
                              &(s->layout.tag), 1)) {
        ESIO_ERROR("Unable to save layout tag", ESIO_ESANITY);
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

// TODO esio_field_open
// TODO use esio_field_{create,open} based upon H5LT_find_dataset

int esio_field_write_double(esio_state s,
                            const char* name,
                            const double *field,
                            int na, int ast, int asz,
                            int nb, int bst, int bsz,
                            int nc, int cst, int csz)
{
    return esio_field_write_internal(s, name, field,
                                     na, ast, asz,
                                     nb, bst, bsz,
                                     nc, cst, csz,
                                     H5T_NATIVE_DOUBLE,
                                     sizeof(double));
}

int esio_field_write_float(esio_state s,
                           const char* name,
                           const float *field,
                           int na, int ast, int asz,
                           int nb, int bst, int bsz,
                           int nc, int cst, int csz)
{
    return esio_field_write_internal(s, name, field,
                                     na, ast, asz,
                                     nb, bst, bsz,
                                     nc, cst, csz,
                                     H5T_NATIVE_FLOAT,
                                     sizeof(float));
}

static
int esio_field_write_internal(esio_state s,
                              const char* name,
                              const void *field,
                              int na, int ast, int asz,
                              int nb, int bst, int bsz,
                              int nc, int cst, int csz,
                              hid_t type_id,
                              size_t type_size)
{
    // Sanity check incoming arguments
    if (s == NULL)        ESIO_ERROR("s == NULL",              ESIO_EINVAL);
    if (s->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EINVAL);
    if (field == NULL)    ESIO_ERROR("field == NULL",          ESIO_EINVAL);
    if (na  < 0)          ESIO_ERROR("na < 0",                 ESIO_EINVAL);
    if (ast < 0)          ESIO_ERROR("ast < 0",                ESIO_EINVAL);
    if (asz < 1)          ESIO_ERROR("asz < 1",                ESIO_EINVAL);
    if (nb  < 0)          ESIO_ERROR("nb < 0",                 ESIO_EINVAL);
    if (bst < 0)          ESIO_ERROR("bst < 0",                ESIO_EINVAL);
    if (bsz < 1)          ESIO_ERROR("bsz < 1",                ESIO_EINVAL);
    if (nc  < 0)          ESIO_ERROR("nc < 0",                 ESIO_EINVAL);
    if (cst < 0)          ESIO_ERROR("cst < 0",                ESIO_EINVAL);
    if (csz < 1)          ESIO_ERROR("csz < 1",                ESIO_EINVAL);

    // TODO Error checking here
    const hid_t dset_id = esio_field_create(s, na, nb, nc, name, type_id);
    (s->layout.field_writer)(dset_id, field,
                             na, ast, asz,
                             nb, bst, bsz,
                             nc, cst, csz,
                             type_id, type_size);
    esio_field_close(dset_id);

    return ESIO_SUCCESS;
}

// **************************************************************
// LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1
// **************************************************************

static
hid_t layout1_filespace_creator(int na, int nb, int nc)
{
    const hsize_t dims[2] = { nb * nc, na };
    return H5Screate_simple(2, dims, NULL);
}

static
hid_t layout1_field_writer(hid_t dset_id, const void *field,
                           int na, int ast, int asz,
                           int nb, int bst, int bsz,
                           int nc, int cst, int csz,
                           hid_t type_id, size_t type_size)
{
    // TODO Error checking here

    (void) nc; // Unused but present for API consistency

    // Create property list for collective write
    const hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // Initialize one-time write details
    const hsize_t stride[2] = { 1, 1 };
    const hsize_t count[2]  = { 1, asz };
    const hid_t   memspace  = H5Screate_simple(2, count, NULL);
    const hid_t   filespace = H5Dget_space(dset_id);

    hsize_t offset[2];
    for (int i = 0; i < csz; ++i)
    {
        for (int j = 0; j < bsz; ++j)
        {
            // Select hyperslab in the file
            offset[0] = (j + bst) + (i + cst) * nb;
            offset[1] = ast;
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
                                offset, stride, count, NULL);

            // Compute in-memory offset to hyperslab's data
            // Note use of type_size when adding to (void *) field
            const size_t moffset = j*asz + i*asz*bsz;
#ifdef __INTEL_COMPILER
#pragma warning(push,disable:1338)
#endif
            const void *p_field = field + (type_size * moffset);
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

            // Write to hyperslab from memory
            const herr_t status = H5Dwrite(dset_id, type_id, memspace,
                                           filespace, plist_id, p_field);
            if (status < 0) {
                H5Sclose(filespace);
                H5Sclose(memspace);
                H5Pclose(plist_id);
                ESIO_ERROR("Write failed", ESIO_EFAILED);
            }
        }
    }

    // Release temporary resources
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);

    return ESIO_SUCCESS;
}
