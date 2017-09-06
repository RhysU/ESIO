//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// ESIO 0.2.0: ExaScale IO library for turbulence simulation restart files
// http://github.com/RhysU/ESIO
//
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
// $Id$

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "esio.h"

#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <mpi.h>

#include "chunksize.h"
#include "error.h"
#include "file-copy.h"
#include "h5utils.h"
#include "layout.h"
#include "metadata.h"
#include "restart-rename.h"
#include "uri.h"
#include "version.h"

//*********************************************************************
// INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL
//*********************************************************************

static
MPI_Comm esio_MPI_Comm_dup_with_name(MPI_Comm comm);

static
hid_t esio_H5P_DATASET_XFER_create(const esio_handle h);

static
hid_t esio_H5P_FILE_ACCESS_create(const esio_handle h);

static
int esio_CONFIGURE_METADATA_CACHING(hid_t plist_id);

static
hid_t esio_field_create(const esio_handle h,
                        const char *name, hid_t type_id,
                        hid_t lcpl_id, hid_t dcpl_id, hid_t dapl_id);

static
hid_t esio_plane_create(const esio_handle h,
                        const char *name, hid_t type_id,
                        hid_t lcpl_id, hid_t dcpl_id, hid_t dapl_id);

static
hid_t esio_line_create(const esio_handle h,
                       const char *name, hid_t type_id,
                       hid_t lcpl_id, hid_t dcpl_id, hid_t dapl_id);

static
int esio_field_close(hid_t dataset_id);

static
int esio_plane_close(hid_t dataset_id);

static
int esio_line_close(hid_t dataset_id);

static
int esio_field_write_internal(const esio_handle h,
                              const char *name,
                              const void *field,
                              int cstride, int bstride, int astride,
                              const char *comment,
                              hid_t type_id);

static
int esio_field_read_internal(const esio_handle h,
                             const char *name,
                             void *field,
                             int cstride, int bstride, int astride,
                             const char *comment,
                             hid_t type_id);

static
int esio_plane_write_internal(const esio_handle h,
                              const char *name,
                              const void *plane,
                              int bstride, int astride,
                              const char *comment,
                              hid_t type_id);

static
int esio_plane_read_internal(const esio_handle h,
                             const char *name,
                             void *plane,
                             int bstride, int astride,
                             const char *comment,
                             hid_t type_id);

static
int esio_line_write_internal(const esio_handle h,
                             const char *name,
                             const void *line,
                             int astride,
                             const char *comment,
                             hid_t type_id);

static
int esio_line_read_internal(const esio_handle h,
                            const char *name,
                            void *line,
                            int astride,
                            const char *comment,
                            hid_t type_id);

// Used in some of the macro-based code generation for conditional arguments
#define WCMTPAR    , const char *comment
#define WCMTARG    , comment
#define RCMTPAR    /* NOP */
#define RCMTARG    , 0

//*********************************************************************
// INTERNAL TYPES INTERNAL TYPES INTERNAL TYPES INTERNAL TYPES INTERNAL
//*********************************************************************

// Bit flags used to control some runtime behavior.
enum {
    FLAG_COLLECTIVE_ENABLED = 1 << 0, //< Should collective IO be used?
    FLAG_CHUNKING_ENABLED   = 1 << 1  //< See features #1246 and #1247
};

struct line_decomp_s {
    int aglobal, astart, alocal;
    int achunk;                   // Cache for when FLAG_CHUNKING_ENABLED
};

struct plane_decomp_s {
    int bglobal, bstart, blocal;
    int aglobal, astart, alocal;
    int bchunk, achunk;           // Cache for when FLAG_CHUNKING_ENABLED
};

struct field_decomp_s {
    int cglobal, cstart, clocal;
    int bglobal, bstart, blocal;
    int aglobal, astart, alocal;
    int cchunk, bchunk, achunk;   // Cache for when FLAG_CHUNKING_ENABLED
};

struct esio_handle_s {
    MPI_Comm  comm;          //< Communicator used for collective calls
    int       comm_rank;     //< Process rank within in MPI communicator
    int       comm_size;     //< Number of ranks within MPI communicator
    MPI_Info  info;          //< Info object used for collective calls
    hid_t     file_id;       //< Active HDF file identifier
    char     *file_path;     //< Active file's canonical path
    int       layout_index;  //< Active field layout_index within HDF5 file
    int       flags;         //< Miscellaneous bit-based flags
    struct line_decomp_s  l; //< Active parallel decomposition for lines
    struct plane_decomp_s p; //< Active parallel decomposition for planes
    struct field_decomp_s f; //< Active parallel decomposition for fields
};

//***************************************************************************
// IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION IMPLEMENTATION
//***************************************************************************

// Lookup table of all the different layout types we understand
static const struct {
    int                      index;
    esio_filespace_creator_t filespace_creator;
    esio_dataset_chunker_t   dataset_chunker;
    esio_field_writer_t      field_writer;
    esio_field_reader_t      field_reader;
} esio_field_layout[] = {
    {
        0,
        &esio_field_layout0_filespace_creator,
        &esio_field_layout0_dataset_chunker,
        &esio_field_layout0_field_writer,
        &esio_field_layout0_field_reader
    },
    {
        1,
        &esio_field_layout1_filespace_creator,
        &esio_field_layout1_dataset_chunker,
        &esio_field_layout1_field_writer,
        &esio_field_layout1_field_reader
    },
    {
        2,
        &esio_field_layout2_filespace_creator,
        &esio_field_layout2_dataset_chunker,
        &esio_field_layout2_field_writer,
        &esio_field_layout2_field_reader
    },
};
static const int esio_field_nlayout = sizeof(esio_field_layout)
                                    / sizeof(esio_field_layout[0]);

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

static
hid_t esio_H5P_DATASET_XFER_create(const esio_handle h)
{
    (void) h; // Unused (for now)

    // Create property list for collective operation
    const hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    if (plist_id < 0) {
        ESIO_ERROR_VAL("Creating plist with type H5P_DATASET_XFER failed",
                       ESIO_ESANITY, -1);
    }

    if (h->flags & FLAG_COLLECTIVE_ENABLED) {
        // Set property list to perform collective operation
        if (H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE) < 0) {
            H5Pclose(plist_id);
            ESIO_ERROR_VAL("Setting IO transfer mode on plist failed",
                           ESIO_ESANITY, -1);
        }
    }

    return plist_id;
}

static
hid_t esio_H5P_FILE_ACCESS_create(const esio_handle h)
{
    // Initialize file access list property identifier
    const hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    if (fapl_id == -1) {
        ESIO_ERROR_VAL("Unable to create fapl_id",
                       ESIO_ESANITY, -1);
    }

    // Set property list collective details
    if (H5Pset_fapl_mpio(fapl_id, h->comm, h->info)) {
        H5Pclose(fapl_id);
        ESIO_ERROR_VAL("Unable to store MPI details in fapl_id",
                       ESIO_ESANITY, -1);
    }

    return fapl_id;
}

static
int esio_CONFIGURE_METADATA_CACHING(hid_t plist_id)
{
    // http://www.hdfgroup.org/pubs/papers/howison_hdf5_lustre_iasds2010.pdf
    // contains a discussion of this logic on pages 3 and 4.  The logic
    // is given (sans error checking) in figure 5.

    // See this question on Hdf-form regarding the H5C_flash_incr__off setting
    // http://mail.hdfgroup.org/pipermail/hdf-forum_hdfgroup.org/2011-February/004199.html

    H5AC_cache_config_t mdc_config;
    mdc_config.version = H5AC__CURR_CACHE_CONFIG_VERSION;

    if (H5Pget_mdc_config(plist_id, &mdc_config) < 0) {
        ESIO_ERROR_VAL("Error calling H5Pget_mdc_config", ESIO_ESANITY, -1);
    }
    mdc_config.evictions_enabled = 0 /* FALSE */;
    mdc_config.incr_mode         = H5C_incr__off;
    mdc_config.flash_incr_mode   = H5C_flash_incr__off;
    mdc_config.decr_mode         = H5C_decr__off;
    if (H5Pset_mdc_config(plist_id, &mdc_config) < 0) {
        ESIO_ERROR_VAL("Error calling H5Pset_mdc_config", ESIO_ESANITY, -1);
    }

    return ESIO_SUCCESS;
}

esio_handle
esio_handle_initialize(MPI_Comm comm)
{
    // Sanity check incoming arguments
    if (comm == MPI_COMM_NULL) {
        ESIO_ERROR_NULL("comm == MPI_COMM_NULL", ESIO_EINVAL);
    }

    // Get number of processors and the local rank within the communicator
    int comm_size;
    ESIO_MPICHKN(MPI_Comm_size(comm, &comm_size));
    int comm_rank;
    ESIO_MPICHKN(MPI_Comm_rank(comm, &comm_rank));

    // Initialize an MPI Info instance
    MPI_Info info;
    ESIO_MPICHKN(MPI_Info_create(&info));

    // Create and initialize ESIO's opaque handler struct
    esio_handle h = calloc(1, sizeof(struct esio_handle_s));
    if (h == NULL) {
        ESIO_ERROR_NULL("failed to allocate space for handle", ESIO_ENOMEM);
    }
    h->comm         = esio_MPI_Comm_dup_with_name(comm);
    h->comm_rank    = comm_rank;
    h->comm_size    = comm_size;
    h->info         = info;
    h->file_id      = -1;
    h->file_path    = NULL;
    h->layout_index = 0;
    h->flags        = FLAG_COLLECTIVE_ENABLED;

    if (h->comm == MPI_COMM_NULL) {
        esio_handle_finalize(h);
        ESIO_ERROR_NULL("Detected MPI_COMM_NULL in h->comm", ESIO_ESANITY);
    }

    return h;
}

esio_handle
esio_handle_initialize_fortran(MPI_Fint fcomm)
{
    // Converting MPI communicators from Fortran to C requires MPI_Comm_f2c
    // See section 16.3.4 of the MPI 2.2 Standard for details
    return esio_handle_initialize(MPI_Comm_f2c(fcomm));
}

int
esio_handle_comm_size(const esio_handle h, int *size)
{
    if (h == NULL)    ESIO_ERROR("h == NULL",    ESIO_EFAULT);
    if (size == NULL) ESIO_ERROR("size == NULL", ESIO_EFAULT);

    *size = h->comm_size;

    return ESIO_SUCCESS;
}

int
esio_handle_comm_rank(const esio_handle h, int *rank)
{
    if (h == NULL)    ESIO_ERROR("h == NULL",    ESIO_EFAULT);
    if (rank == NULL) ESIO_ERROR("rank == NULL", ESIO_EFAULT);

    *rank = h->comm_rank;

    return ESIO_SUCCESS;
}

int
esio_handle_finalize(esio_handle h)
{
    if (h) {
        esio_file_close(h); // Close any open file
        if (h->comm != MPI_COMM_NULL) {
            ESIO_MPICHKR(MPI_Comm_free(&h->comm));
            h->comm = MPI_COMM_NULL;
        }
        if (h->info != MPI_INFO_NULL) {
            ESIO_MPICHKR(MPI_Info_free(&h->info));
            h->info = MPI_INFO_NULL;
        }
        if (h->file_path) {
            free(h->file_path);
            h->file_path = NULL;
        }
        free(h);
    }

    return ESIO_SUCCESS;
}

int esio_field_layout_count()
{
    return esio_field_nlayout;
}

int
esio_field_layout_get(const esio_handle h)
{
    if (h == NULL) {
        // Normal error checking conventions would dictate we invoke
        //   ESIO_ERROR("h == NULL", ESIO_EFAULT);
        // but this routine's return values coincide with the error code range.
        // Instead return the global default layout on NULL input.
        return 0;
    }

    return h->layout_index;
}

int
esio_field_layout_set(esio_handle h, int layout_index)
{
    if (h == NULL) {
        ESIO_ERROR("h == NULL", ESIO_EFAULT);
    }
    if (layout_index < 0) {
        ESIO_ERROR("layout_index < 0", ESIO_EINVAL);
    }
    if (layout_index >= esio_field_nlayout) {
        ESIO_ERROR("layout_index >= esio_field_nlayout", ESIO_EINVAL);
    }

    h->layout_index = layout_index;

    // Changing the layout invalidates the chunksize cache
    h->f.cchunk = h->f.bchunk = h->f.achunk = 0;
    h->p.bchunk = h->p.achunk = 0;
    h->l.achunk = 0;

    return ESIO_SUCCESS;
}

int
esio_file_create(esio_handle h, const char *file, int overwrite)
{
    // Sanity check incoming arguments
    if (h == NULL) {
        ESIO_ERROR("h == NULL", ESIO_EFAULT);
    }
    if (h->file_id != -1) {
        ESIO_ERROR("Cannot create file because previous file not closed",
                   ESIO_EINVAL);
    }
    if (file == NULL) {
        ESIO_ERROR("file == NULL", ESIO_EFAULT);
    }

    // Initialize file creation property list identifier
    const hid_t fcpl_id = H5P_DEFAULT;

    // Initialize file access list property identifier
    const hid_t fapl_id = esio_H5P_FILE_ACCESS_create(h);
    if (fapl_id < 0) {
        ESIO_ERROR("Unable to create fapl_id", ESIO_ESANITY);
    }

    // Set metadata caching options on the file access list property identifier
    if (esio_CONFIGURE_METADATA_CACHING(fapl_id) != ESIO_SUCCESS) {
        H5Pclose(fapl_id);
        ESIO_ERROR("Unable to configure metadata caching", ESIO_ESANITY);
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

    // Clean up temporary HDF5 resources
    H5Pclose(fapl_id);

    // Duplicate the (prefix-less) canonical file name for later use
    // Canonical chosen so that changes in working directory are irrelevant
    h->file_path = canonicalize_file_name(file + scheme_prefix_len(file));
    if (h->file_path == NULL) {
        ESIO_ERROR("failed to allocate space for file_path", ESIO_ENOMEM);
    }

    // File creation successful: update handle
    h->file_id = file_id;

    return ESIO_SUCCESS;
}

int
esio_file_open(esio_handle h, const char *file, int readwrite)
{
    // Sanity check incoming arguments
    if (h == NULL) {
        ESIO_ERROR("h == NULL", ESIO_EFAULT);
    }
    if (h->file_id != -1) {
        ESIO_ERROR("Cannot open new file because previous file not closed",
                   ESIO_EINVAL);
    }
    if (file == NULL) {
        ESIO_ERROR("file == NULL", ESIO_EFAULT);
    }

    // Initialize file access list property identifier
    const hid_t fapl_id = esio_H5P_FILE_ACCESS_create(h);
    if (fapl_id < 0) {
        ESIO_ERROR("Unable to create fapl_id", ESIO_ESANITY);
    }

    // Set metadata caching options on the file access list property identifier
    if (esio_CONFIGURE_METADATA_CACHING(fapl_id) != ESIO_SUCCESS) {
        H5Pclose(fapl_id);
        ESIO_ERROR("Unable to configure metadata caching", ESIO_ESANITY);
    }

    // Initialize access flags
    const unsigned flags = readwrite ? H5F_ACC_RDWR : H5F_ACC_RDONLY;

    // Collectively open the file
    const hid_t file_id = H5Fopen(file, flags, fapl_id);
    if (file_id < 0) {
        H5Pclose(fapl_id);
        ESIO_ERROR("Unable to open existing file", ESIO_EFAILED);
    }

    // Clean up temporary HDF5 resources
    H5Pclose(fapl_id);

    // Duplicate the (prefix-less) canonical file name for later use
    // Canonical chosen so that changes in working directory are irrelevant
    //
    // canonicalize_file_name (aka. realpath) call just below has three issues:
    //   1) A race occurs for file creation versus the canonicalize call
    //   2) Every rank calling canonicalize_file_name hits the file system
    //   3) Memory usage of the returned string isn't promised to be low
    //
    // The messy logic below here is to rectify all three concerns.
    // It is ugly as sin and chatty but pre-allocating everything in advance is
    // and calling realpath(3) is buggy as described in realpath(3)'s man page.

    // Collectively best-effort flush to the operating system
    H5Fflush(file_id, H5F_SCOPE_LOCAL);

    // Last rank canonicalizes the file path, strdups it, and measures length
    // Error information saved but not reported until later.
    const int worker = h->comm_size - 1; // Last rank does work
    int  buf[3] = { /*length*/ -1, /*status*/ ESIO_SUCCESS, /*line*/ -1 };
    char msg[384];
    if (h->comm_rank == worker) {
        static const char msgpat[] = "failed to canonicalize path '%s': %s";
        const char * s = canonicalize_file_name(file + scheme_prefix_len(file));
        if (s == NULL) {
            snprintf(msg, sizeof(msg), msgpat, file, strerror(errno));
            buf[1] = ESIO_EFAILED;
            buf[2] = __LINE__;
        } else if ((h->file_path = strdup(s)) == NULL) {
            snprintf(msg, sizeof(msg), msgpat, file, strerror(ENOMEM));
            buf[1] = ESIO_ENOMEM;
            buf[2] = __LINE__;
        } else if ((buf[0] = (int) strlen(s)) < 1) {
            snprintf(msg, sizeof(msg), msgpat, file, strerror(EINVAL));
            buf[1] = ESIO_ESANITY;
            buf[2] = __LINE__;
        }
        free((void *) s);
    }

    // Broadcast details from rank zero to everyone and return if necessary
    ESIO_MPICHKQ(MPI_Bcast(buf, sizeof(buf)/sizeof(buf[0]), MPI_INT,
                           worker, h->comm));
    if (buf[1] != ESIO_SUCCESS) {
        if (h->comm_rank == worker) {
            esio_error(msg, __FILE__, buf[2], buf[1]);
            free(h->file_path);
        }
        return buf[1];
    }

    // Non-zero ranks allocate space to receive canonical path
    if (h->comm_rank < worker) h->file_path = malloc(buf[0] + 1);
    if (h->file_path == NULL)  buf[1] = ESIO_ENOMEM;

    // Did everybody grab memory?  Otherwise subsequent broadcast crashes...
    ESIO_MPICHKQ(MPI_Allreduce(MPI_IN_PLACE, buf, sizeof(buf)/sizeof(buf[0]),
                               MPI_INT, MPI_MAX, h->comm));
    if (buf[1] != ESIO_SUCCESS) {
        if (h->file_path == NULL) {
            ESIO_ERROR_REPORT("failed to allocate space for canonical path",
                              ESIO_ENOMEM);
        }
        free(h->file_path);
        return buf[1];
    }

    // Broadcast canonical path from rank zero to other ranks
    ESIO_MPICHKQ(MPI_Bcast(h->file_path, buf[0] + 1, MPI_CHAR,
                           worker, h->comm));

    // File creation successful: update handle
    h->file_id = file_id;

    return ESIO_SUCCESS;
}

int esio_file_clone(esio_handle h,
                    const char *srcfile,
                    const char *dstfile,
                    int overwrite)
{
    // Sanity check incoming arguments
    if (h == NULL) {
        ESIO_ERROR("h == NULL", ESIO_EFAULT);
    }
    if (h->file_id != -1) {
        ESIO_ERROR("Cannot create file because previous file not closed",
                   ESIO_EINVAL);
    }
    if (srcfile == NULL) {
        ESIO_ERROR("srcfile == NULL", ESIO_EFAULT);
    }
    if (dstfile == NULL) {
        ESIO_ERROR("dstfile == NULL", ESIO_EFAULT);
    }

    // One rank copies the file synchronously and "broadcasts" the result;
    // The "broadcast" is a summing Allreduce that behaves as useful barrier.
    const int worker = h->comm_size - 1; // Last rank does work
    int status = 0;
    if (h->comm_rank == worker) {
        int prefix_len = scheme_prefix_len(srcfile);
        status = file_copy(srcfile + prefix_len, dstfile + prefix_len,
                           overwrite, 1 /*blockuntilsync*/);
    }
    ESIO_MPICHKQ(MPI_Allreduce(MPI_IN_PLACE, &status, 1/*count*/,
                               MPI_INT, MPI_SUM, h->comm));

    // Bail now if an error occurred
    if (status) return status;

    // All ranks then open the file in readwrite mode
    return esio_file_open(h, dstfile, 1 /* readwrite */);
}

char* esio_file_path(const esio_handle h)
{
    // Sanity check incoming arguments
    if (h == NULL) {
        ESIO_ERROR_NULL("h == NULL", ESIO_EFAULT);
    }
    if (h->file_id == -1) {
        return NULL;
    }
    if (h->file_path == NULL) {
        ESIO_ERROR_NULL("file_id != -1 but file_path == NULL", ESIO_ESANITY);
    }

    // h->file_path has already had any leading URI scheme removed
    char *retval = strdup(h->file_path);
    if (retval == NULL) {
        ESIO_ERROR_NULL("Unable to allocate space for file_path", ESIO_ENOMEM);
    }
    return retval;
}

int esio_file_flush(esio_handle h)
{
    // Sanity check incoming arguments
    if (h == NULL) ESIO_ERROR("h == NULL", ESIO_EFAULT);

    // Flush any currently open file
    if (h->file_id != -1) {
        if (H5Fflush(h->file_id, H5F_SCOPE_GLOBAL) < 0) {
            ESIO_ERROR("Unable to flush file", ESIO_EFAILED);
        }
    }

    return ESIO_SUCCESS;
}

int esio_file_close(esio_handle h)
{
    // Sanity check incoming arguments
    if (h == NULL) ESIO_ERROR("h == NULL", ESIO_EFAULT);

    // Close any currently open file
    if (h->file_id != -1) {

        if (H5Fclose(h->file_id) < 0) {
            ESIO_ERROR("Unable to close file", ESIO_EFAILED);
        }

        if (h->file_path) {
            free(h->file_path);
            h->file_path = NULL;
        }

        // Close successful: update handle
        h->file_id = -1;
    }

    return ESIO_SUCCESS;
}

int esio_file_close_restart(esio_handle h,
                            const char *restart_template,
                            int retain_count)
{
    // Sanity check incoming arguments
    // Template sanity checking is performed within restart_rename as well
    if (h == NULL) {
        ESIO_ERROR("h == NULL", ESIO_EFAULT);
    }
    if (restart_template == NULL) {
        ESIO_ERROR("restart_template == NULL", ESIO_EFAULT);
    }
    if (retain_count < 1) {
        ESIO_ERROR("retain_count < 1", ESIO_EINVAL);
    }
    if (h->file_id == -1) {
        ESIO_ERROR("No file currently open", ESIO_EINVAL);
    }

    // Copy the current file's canonical path and then close the file
    char *src_filename = esio_file_path(h);
    if (src_filename == NULL) {
        ESIO_ERROR("Unable to clone file's canonical path", ESIO_EFAILED);
    }
    const int close_status = esio_file_close(h);
    if (close_status != ESIO_SUCCESS) {
        ESIO_ERROR("Unable to close current restart file", close_status);
    }

    // One rank invokes restart_rename and "broadcasts" the result.
    // The "broadcast" is a summing Allreduce that behaves as useful barrier.
    const int worker = h->comm_size - 1; // Last rank does work
    int status = 0;
    if (h->comm_rank == worker) {
        status = restart_rename(
                src_filename, // No URI scheme prefix munging required
                restart_template + scheme_prefix_len(restart_template),
                retain_count);
    }
    ESIO_MPICHKQ(MPI_Allreduce(MPI_IN_PLACE, &status, 1/*count*/,
                               MPI_INT, MPI_SUM, h->comm));

    // Free temporary memory
    free(src_filename);

    // All ranks return the same status
    if (status != ESIO_SUCCESS) {
        ESIO_ERROR("Failure in esio_file_close_restart", status);
    }
    return ESIO_SUCCESS;
}

static
hid_t esio_field_create(const esio_handle h,
                        const char *name, hid_t type_id,
                        hid_t lcpl_id, hid_t dcpl_id, hid_t dapl_id)
{
    // Sanity check that the handle's layout_index matches our internal table
    if (esio_field_layout[h->layout_index].index != h->layout_index) {
        ESIO_ERROR_VAL("SEVERE: Consistency error in esio_field_layout",
                ESIO_ESANITY, -1);
    }

    // Create the filespace using current layout within handle
    const hid_t filespace
        = (esio_field_layout[h->layout_index].filespace_creator)(h->f.cglobal,
                                                                 h->f.bglobal,
                                                                 h->f.aglobal);
    if (filespace < 0) {
        ESIO_ERROR_VAL("Unable to create filespace", ESIO_ESANITY, -1);
    }

    // Create the dataspace
    const hid_t dset_id = H5Dcreate2(h->file_id, name, type_id, filespace,
                                     lcpl_id, dcpl_id, dapl_id);
    if (dset_id < 0) {
        H5Sclose(filespace);
        ESIO_ERROR_VAL("Unable to create dataspace", ESIO_ESANITY, -1);
    }

    // Stash field's metadata
    if (ESIO_SUCCESS != esio_field_metadata_write(
                h->file_id, name, h->layout_index,
                h->f.cglobal, h->f.bglobal, h->f.aglobal, type_id)) {
        H5Sclose(filespace);
        H5Sclose(dset_id);
        ESIO_ERROR_VAL("Unable to save field metadata", ESIO_EFAILED, -1);
    }

    // Clean up temporary resources
    H5Sclose(filespace);

    return dset_id;
}

static
hid_t esio_plane_create(const esio_handle h,
                        const char *name, hid_t type_id,
                        hid_t lcpl_id, hid_t dcpl_id, hid_t dapl_id)
{
    // Create the filespace
    const hsize_t dims[2] = { h->p.bglobal, h->p.aglobal };
    const hid_t filespace = H5Screate_simple(2, dims, NULL);
    if (filespace < 0) {
        ESIO_ERROR_VAL("Unable to create filespace", ESIO_ESANITY, -1);
    }

    // Create the dataspace
    const hid_t dset_id = H5Dcreate2(h->file_id, name, type_id, filespace,
                                     lcpl_id, dcpl_id, dapl_id);
    if (dset_id < 0) {
        H5Sclose(filespace);
        ESIO_ERROR_VAL("Unable to create dataspace", ESIO_ESANITY, -1);
    }

    // Stash plane's metadata
    if (ESIO_SUCCESS != esio_plane_metadata_write(h->file_id, name,
                                                  h->p.bglobal, h->p.aglobal,
                                                  type_id)) {
        H5Sclose(filespace);
        H5Dclose(dset_id);
        ESIO_ERROR_VAL("Unable to save plane metadata", ESIO_EFAILED, -1);
    }

    // Clean up temporary resources
    H5Sclose(filespace);

    return dset_id;
}

static
hid_t esio_line_create(const esio_handle h,
                       const char *name, hid_t type_id,
                       hid_t lcpl_id, hid_t dcpl_id, hid_t dapl_id)
{
    // Create the filespace
    const hsize_t dims[1] = { h->l.aglobal };
    const hid_t filespace = H5Screate_simple(1, dims, NULL);
    if (filespace < 0) {
        ESIO_ERROR_VAL("Unable to create filespace", ESIO_ESANITY, -1);
    }

    // Create the dataspace
    const hid_t dset_id = H5Dcreate2(h->file_id, name, type_id, filespace,
                                     lcpl_id, dcpl_id, dapl_id);
    if (dset_id < 0) {
        H5Sclose(filespace);
        ESIO_ERROR_VAL("Unable to create dataspace", ESIO_ESANITY, -1);
    }

    // Stash line's metadata
    if (ESIO_SUCCESS != esio_line_metadata_write(h->file_id, name,
                                                 h->l.aglobal, type_id)) {
        H5Sclose(filespace);
        H5Dclose(dset_id);
        ESIO_ERROR_VAL("Unable to save line metadata", ESIO_EFAILED, -1);
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

// *********************************************************************
// PARALLEL DECOMPOSITION DETAILS PARALLEL DECOMPOSITION DETAILS
// *********************************************************************

int
esio_line_establish(esio_handle h,
                    int aglobal, int astart, int alocal)
{
    // Sanity check incoming arguments
    if (h == NULL)   ESIO_ERROR("h == NULL",   ESIO_EFAULT);
    if (aglobal < 1) ESIO_ERROR("aglobal < 1", ESIO_EINVAL);
    if (astart  < 0) ESIO_ERROR("astart < 0",  ESIO_EINVAL);
    if (alocal  < 0) ESIO_ERROR("alocal < 0",  ESIO_EINVAL);

    // Save parallel decomposition in handle
    h->l.aglobal = aglobal;
    h->l.astart  = astart;
    h->l.alocal  = alocal;

    return ESIO_SUCCESS;
}

int
esio_line_established(esio_handle h,
                      int *aglobal, int *astart, int *alocal)
{
    // Sanity check incoming arguments
    if (h == NULL) ESIO_ERROR("h == NULL", ESIO_EFAULT);

    // Retrieve parallel decomposition from handle
    if (aglobal) *aglobal = h->l.aglobal;
    if (astart ) *astart  = h->l.astart;
    if (alocal ) *alocal  = h->l.alocal;

    return ESIO_SUCCESS;
}

int
esio_plane_establish(esio_handle h,
                     int bglobal, int bstart, int blocal,
                     int aglobal, int astart, int alocal)
{
    // Sanity check incoming arguments
    if (h == NULL)   ESIO_ERROR("h == NULL",   ESIO_EFAULT);
    if (bglobal < 1) ESIO_ERROR("bglobal < 1", ESIO_EINVAL);
    if (bstart  < 0) ESIO_ERROR("bstart < 0",  ESIO_EINVAL);
    if (blocal  < 0) ESIO_ERROR("blocal < 0",  ESIO_EINVAL);
    if (aglobal < 1) ESIO_ERROR("aglobal < 1", ESIO_EINVAL);
    if (astart  < 0) ESIO_ERROR("astart < 0",  ESIO_EINVAL);
    if (alocal  < 0) ESIO_ERROR("alocal < 0",  ESIO_EINVAL);

    // Save parallel decomposition in handle
    h->p.bglobal = bglobal;
    h->p.bstart  = bstart;
    h->p.blocal  = blocal;
    h->p.aglobal = aglobal;
    h->p.astart  = astart;
    h->p.alocal  = alocal;

    return ESIO_SUCCESS;
}

int
esio_plane_established(esio_handle h,
                       int *bglobal, int *bstart, int *blocal,
                       int *aglobal, int *astart, int *alocal)
{
    // Sanity check incoming arguments
    if (h == NULL) ESIO_ERROR("h == NULL", ESIO_EFAULT);

    // Retrieve parallel decomposition from handle
    if (bglobal) *bglobal = h->p.bglobal;
    if (bstart ) *bstart  = h->p.bstart;
    if (blocal ) *blocal  = h->p.blocal;
    if (aglobal) *aglobal = h->p.aglobal;
    if (astart ) *astart  = h->p.astart;
    if (alocal ) *alocal  = h->p.alocal;

    return ESIO_SUCCESS;
}

int
esio_field_establish(esio_handle h,
                     int cglobal, int cstart, int clocal,
                     int bglobal, int bstart, int blocal,
                     int aglobal, int astart, int alocal)
{
    // Sanity check incoming arguments
    if (h == NULL)   ESIO_ERROR("h == NULL",   ESIO_EFAULT);
    if (cglobal < 1) ESIO_ERROR("cglobal < 1", ESIO_EINVAL);
    if (cstart  < 0) ESIO_ERROR("cstart < 0",  ESIO_EINVAL);
    if (clocal  < 0) ESIO_ERROR("clocal < 0",  ESIO_EINVAL);
    if (bglobal < 1) ESIO_ERROR("bglobal < 1", ESIO_EINVAL);
    if (bstart  < 0) ESIO_ERROR("bstart < 0",  ESIO_EINVAL);
    if (blocal  < 0) ESIO_ERROR("blocal < 0",  ESIO_EINVAL);
    if (aglobal < 1) ESIO_ERROR("aglobal < 1", ESIO_EINVAL);
    if (astart  < 0) ESIO_ERROR("astart < 0",  ESIO_EINVAL);
    if (alocal  < 0) ESIO_ERROR("alocal < 0",  ESIO_EINVAL);

    // Save parallel decomposition in handle
    h->f.cglobal = cglobal;
    h->f.cstart  = cstart;
    h->f.clocal  = clocal;
    h->f.bglobal = bglobal;
    h->f.bstart  = bstart;
    h->f.blocal  = blocal;
    h->f.aglobal = aglobal;
    h->f.astart  = astart;
    h->f.alocal  = alocal;

    return ESIO_SUCCESS;
}

int
esio_field_established(esio_handle h,
                       int *cglobal, int *cstart, int *clocal,
                       int *bglobal, int *bstart, int *blocal,
                       int *aglobal, int *astart, int *alocal)
{
    // Sanity check incoming arguments
    if (h == NULL) ESIO_ERROR("h == NULL", ESIO_EFAULT);

    // Retrieve parallel decomposition from handle
    if (cglobal) *cglobal = h->f.cglobal;
    if (cstart ) *cstart  = h->f.cstart;
    if (clocal ) *clocal  = h->f.clocal;
    if (bglobal) *bglobal = h->f.bglobal;
    if (bstart ) *bstart  = h->f.bstart;
    if (blocal ) *blocal  = h->f.blocal;
    if (aglobal) *aglobal = h->f.aglobal;
    if (astart ) *astart  = h->f.astart;
    if (alocal ) *alocal  = h->f.alocal;

    return ESIO_SUCCESS;
}

// *********************************************************************
// SIZE SIZEV SIZE SIZEV SIZE SIZEV SIZE SIZEV SIZE SIZEV SIZE SIZEV
// *********************************************************************

int esio_field_size(const esio_handle h,
                    const char *name,
                    int *cglobal, int *bglobal, int *aglobal)
{
    int ncomponents;
    const int status
        = esio_field_sizev(h, name, cglobal, bglobal, aglobal, &ncomponents);
    if (status == ESIO_SUCCESS && ncomponents != 1) {
        ESIO_ERROR("Must retrieve location size using esio_field_sizev",
                   ESIO_EINVAL);
    }
    return status;
}

int esio_field_sizev(const esio_handle h,
                     const char *name,
                     int *cglobal, int *bglobal, int *aglobal,
                     int *ncomponents)
{
    // Sanity check incoming arguments
    if (h == NULL)        ESIO_ERROR("h == NULL",              ESIO_EFAULT);
    if (h->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EFAULT);

    const int status = esio_field_metadata_read(
            h->file_id, name, NULL, cglobal, bglobal, aglobal, ncomponents);
    switch (status) {
        case ESIO_SUCCESS:
        case ESIO_NOTFOUND:  // ESIO_ERROR not called to allow existence query
            return status;
        default:
            ESIO_ERROR("Unable to retrieve field metadata", status);
    }
}

int esio_plane_size(const esio_handle h,
                    const char *name,
                    int *bglobal, int *aglobal)
{
    int ncomponents;
    const int status
        = esio_plane_sizev(h, name, bglobal, aglobal, &ncomponents);
    if (status == ESIO_SUCCESS && ncomponents != 1) {
        ESIO_ERROR("Must retrieve location size using esio_plane_sizev",
                   ESIO_EINVAL);
    }
    return status;
}

int esio_plane_sizev(const esio_handle h,
                     const char *name,
                     int *bglobal, int *aglobal,
                     int *ncomponents)
{
    // Sanity check incoming arguments
    if (h == NULL)        ESIO_ERROR("h == NULL",              ESIO_EFAULT);
    if (h->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EFAULT);

    const int status = esio_plane_metadata_read(
            h->file_id, name, bglobal, aglobal, ncomponents);
    switch (status) {
        case ESIO_SUCCESS:
        case ESIO_NOTFOUND:  // ESIO_ERROR not called to allow existence query
            return status;
        default:
            ESIO_ERROR("Unable to retrieve plane metadata", status);
    }
}

int esio_line_size(const esio_handle h,
                   const char *name,
                   int *aglobal)
{
    int ncomponents;
    const int status = esio_line_sizev(h, name, aglobal, &ncomponents);
    if (status == ESIO_SUCCESS && ncomponents != 1) {
        ESIO_ERROR("Must retrieve location size using esio_line_sizev",
                   ESIO_EINVAL);
    }
    return status;
}

int esio_line_sizev(const esio_handle h,
                    const char *name,
                    int *aglobal,
                    int *ncomponents)
{
    // Sanity check incoming arguments
    if (h == NULL)        ESIO_ERROR("h == NULL",              ESIO_EFAULT);
    if (h->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EFAULT);

    const int status = esio_line_metadata_read(
            h->file_id, name, aglobal, ncomponents);
    switch (status) {
        case ESIO_SUCCESS:
        case ESIO_NOTFOUND:  // ESIO_ERROR not called to allow existence query
            return status;
        default:
            ESIO_ERROR("Unable to retrieve line metadata", status);
    }
}

// *******************************************************************
// FIELD READ WRITE FIELD READ WRITE FIELD READ WRITE FIELD READ WRITE
// *******************************************************************

static
int esio_field_write_internal(const esio_handle h,
                              const char *name,
                              const void *field,
                              int cstride, int bstride, int astride,
                              const char *comment,
                              hid_t type_id)
{
    // Sanity check incoming arguments
    // Strides must be nonnegative because hsize_t is unsigned
    if (h == NULL)        ESIO_ERROR("h == NULL",              ESIO_EFAULT);
    if (h->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EFAULT);
    if (field == NULL && h->f.clocal && h->f.blocal && h->f.alocal) {
                          ESIO_ERROR("field == NULL",           ESIO_EFAULT);
    }
    if (cstride < 0)      ESIO_ERROR("cstride < 0",            ESIO_EINVAL);
    if (bstride < 0)      ESIO_ERROR("bstride < 0",            ESIO_EINVAL);
    if (astride < 0)      ESIO_ERROR("astride < 0",            ESIO_EINVAL);
    // (comment == NULL) is valid input
    if (type_id < 0)      ESIO_ERROR("type_id < 0",            ESIO_EINVAL);
    if (h->f.aglobal == 0)
        ESIO_ERROR("esio_field_establish() never called", ESIO_EINVAL);

    // Provide contiguous defaults whenever the user supplied zero strides.
    // Strides are given in units of type_id; hence astride = 1 is contiguous.
    if (astride == 0) astride = 1;
    if (bstride == 0) bstride = astride * h->f.alocal;
    if (cstride == 0) cstride = bstride * h->f.blocal;

    // Attempt to read metadata for the field (which may or may not exist)
    int layout_index;
    int field_cglobal, field_bglobal, field_aglobal;
    int field_ncomponents;
    const int mstat = esio_field_metadata_read(h->file_id, name,
                                               &layout_index,
                                               &field_cglobal,
                                               &field_bglobal,
                                               &field_aglobal,
                                               &field_ncomponents);

    if (mstat != ESIO_SUCCESS) {
        // Presume field did not exist

        // Determine the chunking parameters to use for dataset creation
        // if they've not already been computed.  Note values are cached!
        if (h->flags & FLAG_CHUNKING_ENABLED && h->f.achunk == 0) {
            const int status = chunksize_field(h->comm, // Expensive
                    h->f.cglobal, h->f.cstart, h->f.clocal, &h->f.cchunk,
                    h->f.bglobal, h->f.bstart, h->f.blocal, &h->f.bchunk,
                    h->f.aglobal, h->f.astart, h->f.alocal, &h->f.achunk);
            if (status != ESIO_SUCCESS) {
                ESIO_ERROR("Error determining chunk size for decomposition",
                        status);
            }
        }

        // Create a dataset creation property list with chunk parameters
        hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        if (dcpl_id < 0) {
            ESIO_ERROR("Error creating dataset creation property list",
                       ESIO_EFAILED);
        }
        if (h->flags & FLAG_CHUNKING_ENABLED) {
            if ((esio_field_layout[h->layout_index].dataset_chunker)(
                        dcpl_id, h->f.cchunk, h->f.bchunk, h->f.achunk) < 0) {
                H5Pclose(dcpl_id);
                ESIO_ERROR("Error setting chunk size information",
                        ESIO_ESANITY);
            }
        }

        // Create dataset and write it with the active field layout
        const hid_t dset_id = esio_field_create(
                h, name, type_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        if (dset_id < 0) {
            ESIO_ERROR("Error creating new field", ESIO_EFAILED);
        }

        // Close creation property list
        H5Pclose(dcpl_id);

        // Obtain appropriate dataset transfer properties
        const hid_t plist_id = esio_H5P_DATASET_XFER_create(h);
        if (plist_id < 0) {
            ESIO_ERROR("Error setting IO transfer properties", ESIO_EFAILED);
        }

        // Write the field using the appropriate layout logic
        const int wstat = (esio_field_layout[h->layout_index].field_writer)(
                plist_id, dset_id, field,
                h->f.cglobal, h->f.cstart, h->f.clocal, cstride,
                h->f.bglobal, h->f.bstart, h->f.blocal, bstride,
                h->f.aglobal, h->f.astart, h->f.alocal, astride,
                type_id);
        if (wstat != ESIO_SUCCESS) {
            esio_field_close(dset_id);
            H5Pclose(plist_id);
            ESIO_ERROR_VAL("Error writing new field", ESIO_EFAILED, wstat);
        }
        H5Pclose(plist_id);

        // Optionally write a comment about the new field
        if (comment && *comment) {
            if (H5Oset_comment(dset_id, comment) < 0) {
                esio_field_close(dset_id);
                ESIO_ERROR("Error setting comment on new field", ESIO_EFAILED);
            }
        }
        esio_field_close(dset_id);

    } else {
        // Field already existed

        // Ensure caller gave correct size information
        if (h->f.cglobal != field_cglobal) {
            ESIO_ERROR("request cglobal mismatch with existing field",
                       ESIO_EINVAL);
        }
        if (h->f.bglobal != field_bglobal) {
            ESIO_ERROR("request bglobal mismatch with existing field",
                       ESIO_EINVAL);
        }
        if (h->f.aglobal != field_aglobal) {
            ESIO_ERROR("request aglobal mismatch with existing field",
                       ESIO_EINVAL);
        }

        // Ensure caller gave type with correct component count
        if (esio_type_ncomponents(type_id) != field_ncomponents) {
            ESIO_ERROR("request ncomponents mismatch with existing field",
                       ESIO_EINVAL);
        }

        // Open the existing field's dataset
        const hid_t dset_id = H5Dopen1(h->file_id, name);
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

        // Obtain appropriate dataset transfer properties
        const hid_t plist_id = esio_H5P_DATASET_XFER_create(h);
        if (plist_id < 0) {
            ESIO_ERROR("Error setting IO transfer properties", ESIO_EFAILED);
        }

        // Overwrite existing data using layout routines per metadata
        const int wstat = (esio_field_layout[layout_index].field_writer)(
                plist_id, dset_id, field,
                h->f.cglobal, h->f.cstart, h->f.clocal, cstride,
                h->f.bglobal, h->f.bstart, h->f.blocal, bstride,
                h->f.aglobal, h->f.astart, h->f.alocal, astride,
                type_id);
        if (wstat != ESIO_SUCCESS) {
            esio_field_close(dset_id);
            H5Pclose(plist_id);
            ESIO_ERROR_VAL("Error overwriting field", ESIO_EFAILED, wstat);
        }
        H5Pclose(plist_id);

        // Optionally write a comment about the field
        if (comment && *comment) {
            if (H5Oset_comment(dset_id, comment) < 0) {
                esio_field_close(dset_id);
                ESIO_ERROR("Error setting comment on existing field",
                        ESIO_EFAILED);
            }
        }
        esio_field_close(dset_id);

    }

    return ESIO_SUCCESS;
}

static
int esio_field_read_internal(const esio_handle h,
                             const char *name,
                             void *field,
                             int cstride, int bstride, int astride,
                             const char *comment,
                             hid_t type_id)
{
    char msg[256]; // message buffer for error handling
    (void) comment; // Present for consistency with esio_field_write_internal
    assert(comment == 0);

    // Sanity check incoming arguments
    // Strides must be nonnegative because hsize_t is unsigned
    if (h == NULL)        ESIO_ERROR("h == NULL",              ESIO_EFAULT);
    if (h->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EFAULT);
    if (field == NULL && h->f.clocal && h->f.blocal && h->f.alocal) {
                          ESIO_ERROR("field == NULL",           ESIO_EFAULT);
    }
    if (cstride < 0)      ESIO_ERROR("cstride < 0",            ESIO_EINVAL);
    if (bstride < 0)      ESIO_ERROR("bstride < 0",            ESIO_EINVAL);
    if (astride < 0)      ESIO_ERROR("astride < 0",            ESIO_EINVAL);
    if (type_id < 0)      ESIO_ERROR("type_id < 0",            ESIO_EINVAL);
    if (h->f.aglobal == 0)
        ESIO_ERROR("esio_field_establish() never called", ESIO_EINVAL);

    // Provide contiguous defaults whenever the user supplied zero strides.
    // Strides are given in units of type_id; hence astride = 1 is contiguous.
    if (astride == 0) astride = 1;
    if (bstride == 0) bstride = astride * h->f.alocal;
    if (cstride == 0) cstride = bstride * h->f.blocal;

    // Read metadata for the field
    int layout_index;
    int field_cglobal, field_bglobal, field_aglobal;
    int field_ncomponents;
    const int status = esio_field_metadata_read(h->file_id, name,
                                                &layout_index,
                                                &field_cglobal,
                                                &field_bglobal,
                                                &field_aglobal,
                                                &field_ncomponents);
    switch (status) {
        case ESIO_SUCCESS:
            break;
        case ESIO_NOTFOUND:
	    snprintf(msg, sizeof(msg), "Field '%s' not found in file", name);
            ESIO_ERROR(msg, status);
        default:
            ESIO_ERROR("Unable to read field's ESIO metadata", status);
    }

    // Ensure caller gave correct size information
    if (h->f.cglobal != field_cglobal) {
        ESIO_ERROR("field read request has incorrect cglobal", ESIO_EINVAL);
    }
    if (h->f.bglobal != field_bglobal) {
        ESIO_ERROR("field read request has incorrect bglobal", ESIO_EINVAL);
    }
    if (h->f.aglobal != field_aglobal) {
        ESIO_ERROR("field read request has incorrect aglobal", ESIO_EINVAL);
    }

    // Ensure caller gave type with correct component count
    if (esio_type_ncomponents(type_id) != field_ncomponents) {
        ESIO_ERROR("request ncomponents mismatch with existing field",
                    ESIO_EINVAL);
    }

    // Open existing dataset
    const hid_t dapl_id = H5P_DEFAULT;
    const hid_t dset_id = H5Dopen2(h->file_id, name, dapl_id);
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

    // Obtain appropriate dataset transfer properties
    const hid_t plist_id = esio_H5P_DATASET_XFER_create(h);
    if (plist_id < 0) {
        H5Dclose(dset_id);
        H5Tclose(field_type_id);
        ESIO_ERROR("Error setting IO transfer properties", ESIO_EFAILED);
    }

    // Read the field based on the metadata's layout_index
    // Note that this means we can read any layout ESIO understands
    // Note that reading does not change the chosen field write layout_index
    (esio_field_layout[layout_index].field_reader)(
            plist_id, dset_id, field,
            h->f.cglobal, h->f.cstart, h->f.clocal, cstride,
            h->f.bglobal, h->f.bstart, h->f.blocal, bstride,
            h->f.aglobal, h->f.astart, h->f.alocal, astride,
            type_id);

    H5Pclose(plist_id);

    // Close dataset
    esio_field_close(dset_id);

    return ESIO_SUCCESS;
}

#define GEN_FIELD_OP(OP,QUAL,TYPE,H5TYPE,CMTPAR,CMTARG)              \
int esio_field_ ## OP ## _ ## TYPE (                                 \
        const esio_handle h,                                         \
        const char *name,                                            \
        QUAL TYPE *field,                                            \
        int cstride, int bstride, int astride                        \
        CMTPAR)                                                      \
{                                                                    \
    return esio_field_ ## OP ## _internal(h, name, field,            \
                                          cstride, bstride, astride  \
                                          CMTARG,                    \
                                          H5TYPE);                   \
}

GEN_FIELD_OP(write, const,       double, H5T_NATIVE_DOUBLE, WCMTPAR, WCMTARG)
GEN_FIELD_OP(read,  /*mutable*/, double, H5T_NATIVE_DOUBLE, RCMTPAR, RCMTARG)

GEN_FIELD_OP(write, const,       float, H5T_NATIVE_FLOAT, WCMTPAR, WCMTARG)
GEN_FIELD_OP(read,  /*mutable*/, float, H5T_NATIVE_FLOAT, RCMTPAR, RCMTARG)

GEN_FIELD_OP(write, const,       int, H5T_NATIVE_INT, WCMTPAR, WCMTARG)
GEN_FIELD_OP(read,  /*mutable*/, int, H5T_NATIVE_INT, RCMTPAR, RCMTARG)


#define GEN_FIELD_OPV(OP,QUAL,TYPE,H5TYPE,CMTPAR,CMTARG)                  \
int esio_field_ ## OP ## v_ ## TYPE(                                      \
        const esio_handle h,                                              \
        const char *name,                                                 \
        QUAL TYPE *field,                                                 \
        int cstride, int bstride, int astride,                            \
        int ncomponents                                                   \
        CMTPAR)                                                           \
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
            h, name, field,                                               \
            (cstride / ncomponents),                                      \
            (bstride / ncomponents),                                      \
            (astride / ncomponents)                                       \
            CMTARG,                                                       \
            array_type_id);                                               \
    H5Tclose(array_type_id);                                              \
    return retval;                                                        \
}

GEN_FIELD_OPV(write, const,       double, H5T_NATIVE_DOUBLE, WCMTPAR, WCMTARG)
GEN_FIELD_OPV(read,  /*mutable*/, double, H5T_NATIVE_DOUBLE, RCMTPAR, RCMTARG)

GEN_FIELD_OPV(write, const,       float, H5T_NATIVE_FLOAT, WCMTPAR, WCMTARG)
GEN_FIELD_OPV(read,  /*mutable*/, float, H5T_NATIVE_FLOAT, RCMTPAR, RCMTARG)

GEN_FIELD_OPV(write, const,       int, H5T_NATIVE_INT, WCMTPAR, WCMTARG)
GEN_FIELD_OPV(read,  /*mutable*/, int, H5T_NATIVE_INT, RCMTPAR, RCMTARG)

// *******************************************************************
// PLANE READ WRITE PLANE READ WRITE PLANE READ WRITE PLANE READ WRITE
// *******************************************************************

static
int esio_plane_write_internal(const esio_handle h,
                              const char *name,
                              const void *plane,
                              int bstride, int astride,
                              const char *comment,
                              hid_t type_id)
{
    // Sanity check incoming arguments
    // Strides must be nonnegative because hsize_t is unsigned
    if (h == NULL)        ESIO_ERROR("h == NULL",              ESIO_EFAULT);
    if (h->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EFAULT);
    if (plane == NULL && h->p.blocal && h->p.alocal) {
                          ESIO_ERROR("plane == NULL",           ESIO_EFAULT);
    }
    if (bstride < 0)      ESIO_ERROR("bstride < 0",            ESIO_EINVAL);
    if (astride < 0)      ESIO_ERROR("astride < 0",            ESIO_EINVAL);
    // (comment == NULL) is valid input
    if (type_id < 0)      ESIO_ERROR("type_id < 0",            ESIO_EINVAL);
    if (h->p.aglobal == 0)
        ESIO_ERROR("esio_plane_establish() never called", ESIO_EINVAL);

    // Provide contiguous defaults whenever the user supplied zero strides.
    // Strides are given in units of type_id; hence astride = 1 is contiguous.
    if (astride == 0) astride = 1;
    if (bstride == 0) bstride = astride * h->p.alocal;

    // Attempt to read metadata for the plane (which may or may not exist)
    int plane_bglobal, plane_aglobal;
    int plane_ncomponents;
    const int mstat = esio_plane_metadata_read(h->file_id, name,
                                               &plane_bglobal,
                                               &plane_aglobal,
                                               &plane_ncomponents);

    hid_t dset_id;
    if (mstat != ESIO_SUCCESS) {
        // Presume plane did not exist

        // Determine the chunking parameters to use for dataset creation
        // if they've not already been computed.  Note values are cached!
        if (h->flags & FLAG_CHUNKING_ENABLED && h->p.achunk == 0) {
            const int status = chunksize_plane(h->comm, // Expensive
                    h->p.bglobal, h->p.bstart, h->p.blocal, &h->p.bchunk,
                    h->p.aglobal, h->p.astart, h->p.alocal, &h->p.achunk);
            if (status != ESIO_SUCCESS) {
                ESIO_ERROR("Error determining chunk size for decomposition",
                        status);
            }
        }

        // Create a dataset creation property list with chunk parameters
        hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        if (dcpl_id < 0) {
            ESIO_ERROR("Error creating dataset creation property list",
                       ESIO_EFAILED);
        }
        if (h->flags & FLAG_CHUNKING_ENABLED) {
            const hsize_t chunksizes[2] = { h->p.bchunk, h->p.achunk };
            if (H5Pset_chunk(dcpl_id, 2, chunksizes) < 0) {
                H5Pclose(dcpl_id);
                ESIO_ERROR("Error setting chunk size information",
                           ESIO_ESANITY);
            }
        }

        // Create the plane
        dset_id = esio_plane_create(
                h, name, type_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        if (dset_id < 0) {
            ESIO_ERROR("Error creating new plane", ESIO_EFAILED);
        }

        // Close creation property list
        H5Pclose(dcpl_id);

    } else {
        // Plane already existed

        // Ensure caller gave correct size information
        if (h->p.bglobal != plane_bglobal) {
            ESIO_ERROR("request bglobal mismatch with existing plane",
                       ESIO_EINVAL);
        }
        if (h->p.aglobal != plane_aglobal) {
            ESIO_ERROR("request aglobal mismatch with existing plane",
                       ESIO_EINVAL);
        }

        // Ensure caller gave type with correct component count
        if (esio_type_ncomponents(type_id) != plane_ncomponents) {
            ESIO_ERROR("request ncomponents mismatch with existing plane",
                       ESIO_EINVAL);
        }

        // Open the existing plane's dataset
        dset_id = H5Dopen1(h->file_id, name);
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

    // Obtain appropriate dataset transfer properties
    const hid_t plist_id = esio_H5P_DATASET_XFER_create(h);
    if (plist_id < 0) {
        esio_plane_close(dset_id);
        ESIO_ERROR("Error setting IO transfer properties", ESIO_EFAILED);
    }

    // Write plane
    const int wstat = esio_plane_writer(
            plist_id, dset_id, plane,
            h->p.bglobal, h->p.bstart, h->p.blocal, bstride,
            h->p.aglobal, h->p.astart, h->p.alocal, astride,
            type_id);
    if (wstat != ESIO_SUCCESS) {
        esio_plane_close(dset_id);
        H5Pclose(plist_id);
        ESIO_ERROR_VAL("Error writing plane", ESIO_EFAILED, wstat);
    }
    H5Pclose(plist_id);

    // Optionally write a comment about the plane
    if (comment && *comment) {
        if (H5Oset_comment(dset_id, comment) < 0) {
            esio_plane_close(dset_id);
            ESIO_ERROR("Error setting comment on plane", ESIO_EFAILED);
        }
    }
    esio_plane_close(dset_id);

    return ESIO_SUCCESS;
}

static
int esio_plane_read_internal(const esio_handle h,
                             const char *name,
                             void *plane,
                             int bstride, int astride,
                             const char *comment,
                             hid_t type_id)
{
    char msg[256]; // message error buffer
    (void) comment; // Present for consistency with esio_field_write_internal
    assert(comment == 0);

    // Sanity check incoming arguments
    // Strides must be nonnegative because hsize_t is unsigned
    if (h == NULL)        ESIO_ERROR("h == NULL",              ESIO_EFAULT);
    if (h->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EFAULT);
    if (plane == NULL && h->p.blocal && h->p.alocal) {
                          ESIO_ERROR("plane == NULL",           ESIO_EFAULT);
    }
    if (bstride < 0)      ESIO_ERROR("bstride < 0",            ESIO_EINVAL);
    if (astride < 0)      ESIO_ERROR("astride < 0",            ESIO_EINVAL);
    if (type_id < 0)      ESIO_ERROR("type_id < 0",            ESIO_EINVAL);
    if (h->p.aglobal == 0)
        ESIO_ERROR("esio_plane_establish() never called", ESIO_EINVAL);

    // Provide contiguous defaults whenever the user supplied zero strides.
    // Strides are given in units of type_id; hence astride = 1 is contiguous.
    if (astride == 0) astride = 1;
    if (bstride == 0) bstride = astride * h->p.alocal;

    // Read metadata for the plane
    int plane_bglobal, plane_aglobal;
    int plane_ncomponents;
    const int status = esio_plane_metadata_read(h->file_id, name,
                                                &plane_bglobal,
                                                &plane_aglobal,
                                                &plane_ncomponents);
    switch (status) {
        case ESIO_SUCCESS:
            break;
        case ESIO_NOTFOUND:
	    snprintf(msg, sizeof(msg), "Plane '%s' not found in file", name);
            ESIO_ERROR(msg, status);
        default:
            ESIO_ERROR("Unable to read plane's ESIO metadata", status);
    }

    // Ensure caller gave correct size information
    if (h->p.bglobal != plane_bglobal) {
        ESIO_ERROR("plane read request has incorrect bglobal", ESIO_EINVAL);
    }
    if (h->p.aglobal != plane_aglobal) {
        ESIO_ERROR("plane read request has incorrect aglobal", ESIO_EINVAL);
    }

    // Ensure caller gave type with correct component count
    if (esio_type_ncomponents(type_id) != plane_ncomponents) {
        ESIO_ERROR("request ncomponents mismatch with existing plane",
                    ESIO_EINVAL);
    }

    // Obtain appropriate dataset transfer properties
    const hid_t plist_id = esio_H5P_DATASET_XFER_create(h);
    if (plist_id < 0) {
        ESIO_ERROR("Error setting IO transfer properties", ESIO_EFAILED);
    }

    // Open existing dataset
    const hid_t dapl_id = H5P_DEFAULT;
    const hid_t dset_id = H5Dopen2(h->file_id, name, dapl_id);
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
        H5Pclose(plist_id);
        ESIO_ERROR("request type not convertible to existing plane type",
                    ESIO_EINVAL);
    }
    H5Tclose(plane_type_id);

    // Read plane
    const int rstat = esio_plane_reader(
            plist_id, dset_id, plane,
            h->p.bglobal, h->p.bstart, h->p.blocal, bstride,
            h->p.aglobal, h->p.astart, h->p.alocal, astride,
            type_id);
    if (rstat != ESIO_SUCCESS) {
        esio_plane_close(dset_id);
        H5Pclose(plist_id);
        ESIO_ERROR_VAL("Error reading plane", ESIO_EFAILED, rstat);
    }
    esio_plane_close(dset_id);
    H5Pclose(plist_id);

    return ESIO_SUCCESS;
}

#define GEN_PLANE_OP(OP,QUAL,TYPE,H5TYPE,CMTPAR,CMTARG)     \
int esio_plane_ ## OP ## _ ## TYPE (                        \
        const esio_handle h,                                \
        const char *name,                                   \
        QUAL TYPE *plane,                                   \
        int bstride, int astride                            \
        CMTPAR)                                             \
{                                                           \
    return esio_plane_ ## OP ## _internal(h, name, plane,   \
                                          bstride, astride  \
                                          CMTARG,           \
                                          H5TYPE);          \
}

GEN_PLANE_OP(write, const,       double, H5T_NATIVE_DOUBLE, WCMTPAR, WCMTARG)
GEN_PLANE_OP(read,  /*mutable*/, double, H5T_NATIVE_DOUBLE, RCMTPAR, RCMTARG)

GEN_PLANE_OP(write, const,       float, H5T_NATIVE_FLOAT, WCMTPAR, WCMTARG)
GEN_PLANE_OP(read,  /*mutable*/, float, H5T_NATIVE_FLOAT, RCMTPAR, RCMTARG)

GEN_PLANE_OP(write, const,       int, H5T_NATIVE_INT, WCMTPAR, WCMTARG)
GEN_PLANE_OP(read,  /*mutable*/, int, H5T_NATIVE_INT, RCMTPAR, RCMTARG)

#define GEN_PLANE_OPV(OP,QUAL,TYPE,H5TYPE,CMTPAR,CMTARG)                  \
int esio_plane_ ## OP ## v_ ## TYPE(                                      \
        const esio_handle h,                                              \
        const char *name,                                                 \
        QUAL TYPE *plane,                                                 \
        int bstride, int astride,                                         \
        int ncomponents                                                   \
        CMTPAR)                                                           \
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
            h, name, plane,                                               \
            (bstride / ncomponents),                                      \
            (astride / ncomponents)                                       \
            CMTARG,                                                       \
            array_type_id);                                               \
    H5Tclose(array_type_id);                                              \
    return retval;                                                        \
}

GEN_PLANE_OPV(write, const,       double, H5T_NATIVE_DOUBLE, WCMTPAR, WCMTARG)
GEN_PLANE_OPV(read,  /*mutable*/, double, H5T_NATIVE_DOUBLE, RCMTPAR, RCMTARG)

GEN_PLANE_OPV(write, const,       float, H5T_NATIVE_FLOAT, WCMTPAR, WCMTARG)
GEN_PLANE_OPV(read,  /*mutable*/, float, H5T_NATIVE_FLOAT, RCMTPAR, RCMTARG)

GEN_PLANE_OPV(write, const,       int, H5T_NATIVE_INT, WCMTPAR, WCMTARG)
GEN_PLANE_OPV(read,  /*mutable*/, int, H5T_NATIVE_INT, RCMTPAR, RCMTARG)

// *******************************************************************
// LINE READ WRITE LINE READ WRITE LINE READ WRITE LINE READ WRITE
// *******************************************************************

static
int esio_line_write_internal(const esio_handle h,
                             const char *name,
                             const void *line,
                             int astride,
                             const char *comment,
                             hid_t type_id)
{
    // Sanity check incoming arguments
    // Strides must be nonnegative because hsize_t is unsigned
    if (h == NULL)        ESIO_ERROR("h == NULL",              ESIO_EFAULT);
    if (h->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EFAULT);
    if (line == NULL && h->l.alocal) {
                          ESIO_ERROR("line == NULL",           ESIO_EFAULT);
    }
    if (astride < 0)      ESIO_ERROR("astride < 0",            ESIO_EINVAL);
    if (type_id < 0)      ESIO_ERROR("type_id < 0",            ESIO_EINVAL);
    if (h->l.aglobal == 0)
        ESIO_ERROR("esio_line_establish() never called", ESIO_EINVAL);

    // Provide contiguous defaults whenever the user supplied zero strides.
    // Strides are given in units of type_id; hence astride = 1 is contiguous.
    if (astride == 0) astride = 1;

    // Attempt to read metadata for the line (which may or may not exist)
    int line_aglobal, line_ncomponents;
    const int mstat = esio_line_metadata_read(h->file_id, name,
                                              &line_aglobal,
                                              &line_ncomponents);

    hid_t dset_id;
    if (mstat != ESIO_SUCCESS) {
        // Presume line did not already exist

        // Determine the chunking parameters to use for dataset creation
        // if they've not already been computed.  Note values are cached!
        if (h->flags & FLAG_CHUNKING_ENABLED && h->l.achunk == 0) {
            const int status = chunksize_line(h->comm, // Expensive
                    h->l.aglobal, h->l.astart, h->l.alocal, &h->l.achunk);
            if (status != ESIO_SUCCESS) {
                ESIO_ERROR("Error determining chunk size for decomposition",
                           status);
            }
        }

        // Create a dataset creation property list with chunk parameters
        hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        if (dcpl_id < 0) {
            ESIO_ERROR("Error creating dataset creation property list",
                       ESIO_EFAILED);
        }
        if (h->flags & FLAG_CHUNKING_ENABLED) {
            const hsize_t chunksizes[1] = { h->l.achunk };
            if (H5Pset_chunk(dcpl_id, 1, chunksizes) < 0) {
                H5Pclose(dcpl_id);
                ESIO_ERROR("Error setting chunk size information",
                           ESIO_ESANITY);
            }
        }

        // Create the line
        dset_id = esio_line_create(
                h, name, type_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        if (dset_id < 0) {
            ESIO_ERROR("Error creating new line", ESIO_EFAILED);
        }

        // Close creation property list
        H5Pclose(dcpl_id);

    } else {
        // Line already existed

        // Ensure caller gave correct size information
        if (h->l.aglobal != line_aglobal) {
            ESIO_ERROR("request aglobal mismatch with existing line",
                       ESIO_EINVAL);
        }

        // Ensure caller gave type with correct component count
        if (esio_type_ncomponents(type_id) != line_ncomponents) {
            ESIO_ERROR("request ncomponents mismatch with existing line",
                       ESIO_EINVAL);
        }

        // Open the existing line's dataset
        dset_id = H5Dopen1(h->file_id, name);
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

    // Obtain appropriate dataset transfer properties
    const hid_t plist_id = esio_H5P_DATASET_XFER_create(h);
    if (plist_id < 0) {
        ESIO_ERROR("Error setting IO transfer properties", ESIO_EFAILED);
    }

    // Write line
    const int wstat = esio_line_writer(
            plist_id, dset_id, line,
            h->l.aglobal, h->l.astart, h->l.alocal, astride, type_id);
    if (wstat != ESIO_SUCCESS) {
        esio_line_close(dset_id);
        H5Pclose(plist_id);
        ESIO_ERROR_VAL("Error writing line", ESIO_EFAILED, wstat);
    }
    H5Pclose(plist_id);

    // Optionally write a comment about the line
    if (comment && *comment) {
        if (H5Oset_comment(dset_id, comment) < 0) {
            esio_line_close(dset_id);
            ESIO_ERROR("Error setting comment on line", ESIO_EFAILED);
        }
    }
    esio_line_close(dset_id);

    return ESIO_SUCCESS;
}

static
int esio_line_read_internal(const esio_handle h,
                            const char *name,
                            void *line,
                            int astride,
                            const char *comment,
                            hid_t type_id)
{
    char msg[256];   // message buffer for error handling
    (void) comment; // Present for consistency with esio_field_write_internal
    assert(comment == 0);

    // Sanity check incoming arguments
    // Strides must be nonnegative because hsize_t is unsigned
    if (h == NULL)        ESIO_ERROR("h == NULL",              ESIO_EFAULT);
    if (h->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EFAULT);
    if (line == NULL && h->l.alocal) {
                          ESIO_ERROR("line == NULL",           ESIO_EFAULT);
    }
    if (astride < 0)      ESIO_ERROR("astride < 0",            ESIO_EINVAL);
    if (type_id < 0)      ESIO_ERROR("type_id < 0",            ESIO_EINVAL);
    if (h->l.aglobal == 0) {
        ESIO_ERROR("esio_line_establish() never called", ESIO_EINVAL);
    }

    // Provide contiguous defaults whenever the user supplied zero strides.
    // Strides are given in units of type_id; hence astride = 1 is contiguous.
    if (astride == 0) astride = 1;

    // Read metadata for the line
    int line_aglobal, line_ncomponents;
    const int status = esio_line_metadata_read(h->file_id, name,
                                               &line_aglobal,
                                               &line_ncomponents);
    switch (status) {
        case ESIO_SUCCESS:
            break;
        case ESIO_NOTFOUND:
	    snprintf(msg, sizeof(msg), "Line '%s' not found in file", name);
            ESIO_ERROR(msg, status);
        default:
            ESIO_ERROR("Unable to read line's ESIO metadata", status);
    }

    // Ensure caller gave correct size information
    if (h->l.aglobal != line_aglobal) {
        ESIO_ERROR("line read request has incorrect aglobal", ESIO_EINVAL);
    }

    // Ensure caller gave type with correct component count
    if (esio_type_ncomponents(type_id) != line_ncomponents) {
        ESIO_ERROR("request ncomponents mismatch with existing line",
                    ESIO_EINVAL);
    }

    // Open existing dataset
    const hid_t dapl_id = H5P_DEFAULT;
    const hid_t dset_id = H5Dopen2(h->file_id, name, dapl_id);
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

    // Obtain appropriate dataset transfer properties
    const hid_t plist_id = esio_H5P_DATASET_XFER_create(h);
    if (plist_id < 0) {
        ESIO_ERROR("Error setting IO transfer properties", ESIO_EFAILED);
    }

    // Read line
    const int rstat = esio_line_reader(
            plist_id, dset_id, line,
            h->l.aglobal, h->l.astart, h->l.alocal, astride, type_id);
    if (rstat != ESIO_SUCCESS) {
        esio_line_close(dset_id);
        H5Pclose(plist_id);
        ESIO_ERROR_VAL("Error reading line", ESIO_EFAILED, rstat);
    }
    esio_line_close(dset_id);
    H5Pclose(plist_id);

    return ESIO_SUCCESS;
}

#define GEN_LINE_OP(OP,QUAL,TYPE,H5TYPE,CMTPAR,CMTARG)  \
int esio_line_ ## OP ## _ ## TYPE (                     \
        const esio_handle h,                            \
        const char *name,                               \
        QUAL TYPE *line,                                \
        int astride                                     \
        CMTPAR)                                         \
{                                                       \
    return esio_line_ ## OP ## _internal(h, name, line, \
                                         astride        \
                                         CMTARG,        \
                                         H5TYPE);       \
}

GEN_LINE_OP(write, const,       double, H5T_NATIVE_DOUBLE, WCMTPAR, WCMTARG)
GEN_LINE_OP(read,  /*mutable*/, double, H5T_NATIVE_DOUBLE, RCMTPAR, RCMTARG)

GEN_LINE_OP(write, const,       float, H5T_NATIVE_FLOAT, WCMTPAR, WCMTARG)
GEN_LINE_OP(read,  /*mutable*/, float, H5T_NATIVE_FLOAT, RCMTPAR, RCMTARG)

GEN_LINE_OP(write, const,       int, H5T_NATIVE_INT, WCMTPAR, WCMTARG)
GEN_LINE_OP(read,  /*mutable*/, int, H5T_NATIVE_INT, RCMTPAR, RCMTARG)

#define GEN_LINE_OPV(OP,QUAL,TYPE,H5TYPE,CMTPAR,CMTARG)                   \
int esio_line_ ## OP ## v_ ## TYPE(                                       \
        const esio_handle h,                                              \
        const char *name,                                                 \
        QUAL TYPE *line,                                                  \
        int astride,                                                      \
        int ncomponents                                                   \
        CMTPAR)                                                           \
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
            h, name, line,                                                \
            (astride / ncomponents)                                       \
            CMTARG,                                                       \
            array_type_id);                                               \
    H5Tclose(array_type_id);                                              \
    return retval;                                                        \
}

GEN_LINE_OPV(write, const,       double, H5T_NATIVE_DOUBLE, WCMTPAR, WCMTARG)
GEN_LINE_OPV(read,  /*mutable*/, double, H5T_NATIVE_DOUBLE, RCMTPAR, RCMTARG)

GEN_LINE_OPV(write, const,       float, H5T_NATIVE_FLOAT, WCMTPAR, WCMTARG)
GEN_LINE_OPV(read,  /*mutable*/, float, H5T_NATIVE_FLOAT, RCMTPAR, RCMTARG)

GEN_LINE_OPV(write, const,       int, H5T_NATIVE_INT, WCMTPAR, WCMTARG)
GEN_LINE_OPV(read,  /*mutable*/, int, H5T_NATIVE_INT, RCMTPAR, RCMTARG)

// *********************************************************************
// ATTRIBUTE ATTRIBUTE ATTRIBUTE ATTRIBUTE ATTRIBUTE ATTRIBUTE ATTRIBUTE
// *********************************************************************

int esio_attribute_sizev(const esio_handle h,
                         const char *location,
                         const char *name,
                         int *ncomponents)
{
    // Sanity check incoming arguments
    if (h == NULL)        ESIO_ERROR("h == NULL",              ESIO_EFAULT);
    if (h->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (location == NULL) ESIO_ERROR("location == NULL",           ESIO_EFAULT);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EFAULT);

    // Attempt to retrieve information on the attribute
    int rank;
    hsize_t dims[H5S_MAX_RANK]; // Oversized to protect smashing the stack
    H5T_class_t type_class;
    size_t type_size;
    DISABLE_HDF5_ERROR_HANDLER(one)
    const herr_t err = esio_H5LTget_attribute_ndims_info(h->file_id,
                                                         location,
                                                         name,
                                                         &rank, dims,
                                                         &type_class,
                                                         &type_size);
    ENABLE_HDF5_ERROR_HANDLER(one)
    if (err == H5E_NOTFOUND) {
        return ESIO_NOTFOUND;  // ESIO_ERROR not called to allow existence query
    } else if (err < 0) {
        ESIO_ERROR("Failure querying attribute at location", ESIO_EFAILED);
    }
    if (rank != 1) {
        ESIO_ERROR("Attribute rank != 1 unsupported", ESIO_EFAILED);
    }

    // Ensure we won't overflow the return value type
    if (dims[0] > INT_MAX) {
        ESIO_ERROR("Attribute size > INT_MAX", ESIO_ESANITY);
    }

    // Return requested information to the caller
    if (ncomponents) *ncomponents = dims[0];

    return ESIO_SUCCESS;
}

#define GEN_ATTRIBUTE_WRITEV(TYPE)                                            \
int esio_attribute_writev_##TYPE(const esio_handle h,                         \
                                 const char *location,                        \
                                 const char *name,                            \
                                 const TYPE *value,                           \
                                 int ncomponents)                             \
{                                                                             \
    /* Sanity check incoming arguments */                                     \
    if (h == NULL)        ESIO_ERROR("h == NULL",              ESIO_EFAULT);  \
    if (h->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);  \
    if (location == NULL) ESIO_ERROR("location == NULL",       ESIO_EFAULT);  \
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EFAULT);  \
    if (value == NULL)    ESIO_ERROR("value == NULL",          ESIO_EFAULT);  \
    if (ncomponents < 1)  ESIO_ERROR("ncomponents < 1",        ESIO_EINVAL);  \
                                                                              \
    const herr_t err = H5LTset_attribute_##TYPE(                              \
            h->file_id, location, name, value, ncomponents);                  \
    if (err < 0) {                                                            \
        ESIO_ERROR("unable to write attribute at location", ESIO_EFAILED);    \
    }                                                                         \
                                                                              \
    return ESIO_SUCCESS;                                                      \
}

#define GEN_ATTRIBUTE_READV(TYPE)                                             \
int esio_attribute_readv_##TYPE(const esio_handle h,                          \
                                const char *location,                         \
                                const char *name,                             \
                                TYPE *value,                                  \
                                int ncomponents)                              \
{                                                                             \
    char msg[256]; /* message buffer for errors */                            \
                                                                              \
    /* Sanity check incoming arguments */                                     \
    if (h == NULL)        ESIO_ERROR("h == NULL",              ESIO_EFAULT);  \
    if (h->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);  \
    if (location == NULL) ESIO_ERROR("location == NULL",       ESIO_EFAULT);  \
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EFAULT);  \
    if (value == NULL)    ESIO_ERROR("value == NULL",          ESIO_EFAULT);  \
    if (ncomponents < 1)  ESIO_ERROR("ncomponents < 1",        ESIO_EINVAL);  \
                                                                              \
    /* Attempt to retrieve information on the attribute */                    \
    int rank;                                                                 \
    hsize_t dims[H5S_MAX_RANK]; /* Oversized protects the stack */            \
    H5T_class_t type_class;                                                   \
    size_t type_size;                                                         \
    DISABLE_HDF5_ERROR_HANDLER(one)                                           \
    const herr_t err1 = esio_H5LTget_attribute_ndims_info(                    \
            h->file_id, location, name, &rank, dims,                          \
            &type_class, &type_size);                                         \
    ENABLE_HDF5_ERROR_HANDLER(one)                                            \
    if (err1 == H5E_NOTFOUND) {                                               \
        snprintf(msg, sizeof(msg), "Attribute '%s' not found at location",    \
                 name);                                                       \
        ESIO_ERROR(msg, ESIO_EINVAL);                                         \
    } else if (err1 < 0) {                                                    \
        ESIO_ERROR("unable to interrogate attribute at location",             \
                    ESIO_EINVAL);                                             \
    }                                                                         \
                                                                              \
    /* Ensure the user is getting what he/she requested */                    \
    if ((hsize_t) ncomponents != dims[0]) {                                   \
        ESIO_ERROR("requested ncomponents mismatch with stored ncomponents",  \
                   ESIO_EINVAL);                                              \
    }                                                                         \
                                                                              \
    /* Read the attribute's data */                                           \
    const herr_t err2 = H5LTget_attribute_##TYPE(                             \
            h->file_id, location, name, value);                               \
    if (err2 < 0) {                                                           \
        ESIO_ERROR("unable to retrieve attribute at location",                \
                ESIO_EFAILED);                                                \
    }                                                                         \
                                                                              \
    return ESIO_SUCCESS;                                                      \
}

#define GEN_ATTRIBUTE_WRITE(TYPE)                                             \
int esio_attribute_write_##TYPE(const esio_handle h,                          \
                                const char *location,                         \
                                const char *name,                             \
                                const TYPE *value)                            \
{                                                                             \
    return esio_attribute_writev_##TYPE(h, location, name, value, 1);         \
}

#define GEN_ATTRIBUTE_READ(TYPE)                                              \
int esio_attribute_read_##TYPE(const esio_handle h,                           \
                               const char *location,                          \
                               const char *name,                              \
                               TYPE *value)                                   \
{                                                                             \
    return esio_attribute_readv_##TYPE(h, location, name, value, 1);          \
}

GEN_ATTRIBUTE_WRITEV(double)
GEN_ATTRIBUTE_WRITE(double)
GEN_ATTRIBUTE_READV(double)
GEN_ATTRIBUTE_READ(double)

GEN_ATTRIBUTE_WRITEV(float)
GEN_ATTRIBUTE_WRITE(float)
GEN_ATTRIBUTE_READV(float)
GEN_ATTRIBUTE_READ(float)

GEN_ATTRIBUTE_WRITEV(int)
GEN_ATTRIBUTE_WRITE(int)
GEN_ATTRIBUTE_READV(int)
GEN_ATTRIBUTE_READ(int)

// *********************************************************************
// STRING STRING STRING STRING STRING STRING STRING STRING STRING STRING
// *********************************************************************

int esio_string_set(const esio_handle h,
                    const char *location,
                    const char *name,
                    const char *value)
{
    // Sanity check incoming arguments
    if (h == NULL)        ESIO_ERROR("h == NULL",              ESIO_EFAULT);
    if (h->file_id == -1) ESIO_ERROR("No file currently open", ESIO_EINVAL);
    if (location == NULL) ESIO_ERROR("location == NULL",       ESIO_EFAULT);
    if (name == NULL)     ESIO_ERROR("name == NULL",           ESIO_EFAULT);
    if (value == NULL)    ESIO_ERROR("value == NULL",          ESIO_EFAULT);

    const herr_t err = H5LTset_attribute_string(
            h->file_id, location, name, value);
    if (err < 0) {
        ESIO_ERROR("unable to write attribute at location", ESIO_EFAILED);
    }

    return ESIO_SUCCESS;
}

char* esio_string_get(const esio_handle h,
                      const char *location,
                      const char *name)
{
    // Sanity check incoming arguments
    if (h == NULL)
        ESIO_ERROR_NULL("h == NULL",              ESIO_EFAULT);
    if (h->file_id == -1)
        ESIO_ERROR_NULL("No file currently open", ESIO_EINVAL);
    if (location == NULL)
        ESIO_ERROR_NULL("location == NULL",       ESIO_EFAULT);
    if (name == NULL)
        ESIO_ERROR_NULL("name == NULL",           ESIO_EFAULT);

    // Silence HDF5 errors while querying the attribute
    DISABLE_HDF5_ERROR_HANDLER(one)

    // Attempt to open the requested attribute
    const hid_t aid = H5Aopen_by_name(
            h->file_id, location, name, H5P_DEFAULT, H5P_DEFAULT);
    if (aid < 0) {
        ENABLE_HDF5_ERROR_HANDLER(one)
        if (0 == H5Aexists_by_name(h->file_id, location, name, H5P_DEFAULT)) {
            return NULL;
        } else {
            ESIO_ERROR_NULL("unable to interrogate attribute at location",
                            ESIO_EINVAL);
        }
    }

    // Open type associated with the attribute
    const hid_t tid = H5Aget_type(aid);
    if (tid < 0) {
        ENABLE_HDF5_ERROR_HANDLER(one)
        H5Aclose(aid);
        ESIO_ERROR_NULL("unable to interrogate attribute type", ESIO_EINVAL);
    }

    // Check that we're dealing with a string
    if (H5T_STRING != H5Tget_class(tid)) {
        ENABLE_HDF5_ERROR_HANDLER(one)
        H5Tclose(tid);
        H5Aclose(aid);
        ESIO_ERROR_NULL("attribute is not a string", ESIO_EINVAL);
    }

    // Below here HDF5 calls should be succeeding
    ENABLE_HDF5_ERROR_HANDLER(one)

    // Retrieve the string length
    const size_t size = H5Tget_size(tid);
    if (size == 0) {
        H5Tclose(tid);
        H5Aclose(aid);
        ESIO_ERROR_NULL("unable to obtain string length",
                ESIO_EINVAL);
    }

    // Allocate storage for the string (already includes null termination)
    char * retval = malloc(size);
    if (retval == NULL) {
        H5Tclose(tid);
        H5Aclose(aid);
        ESIO_ERROR_NULL("Unable to allocate storage for string", ESIO_ENOMEM);
    }

    // Read the string's data
    const herr_t err = H5Aread(aid, tid, retval);
    if (err < 0) {
        H5Tclose(tid);
        H5Aclose(aid);
        ESIO_ERROR_NULL("unable to retrieve string", ESIO_EFAILED);
    }

    // Close the attribute type and the attribute
    H5Tclose(tid);
    H5Aclose(aid);

    // Return the newly allocated string.
    // The caller MUST free the memory.
    return retval;
}
