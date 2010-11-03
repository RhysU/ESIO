//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// esio 0.0.1: ExaScale IO library for turbulence simulation restart files
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

#ifndef __ESIO_ESIO_H
#define __ESIO_ESIO_H

#include <mpi.h>

/** @file
 * Provides ESIO's C-based public API following the library's
 * \ref concepts "usage concepts".  All methods in this file
 * invoke ESIO's \ref errorC "error handling" mechanisms on failure.
 */

#ifdef __cplusplus
extern "C" {
#endif

/** An opaque type following ESIO's \ref conceptshandles "handle concept". */
typedef struct esio_state_s *esio_state;

/**
 * \name Initializing and finalizing a state handle
 * See \ref conceptshandles "handles" for the associated semantics.
 */
/*\@{*/

/**
 * Initialize a handle against the given MPI communicator.
 * The handle must be finalized using esio_finalize() to avoid
 * resource leaks.
 *
 * \param comm MPI communicator (e.g. \c MPI_COMM_WORLD) used to determine
 *             the parallel scope of the handle.  All ESIO calls employing
 *             the handle must be made collectively on \c comm.
 * \return A new handle on success.  Otherwise \c NULL.
 */
esio_state esio_initialize(MPI_Comm comm);

/**
 * Finalize a handle.
 * Finalizing a handle automatically closes any associated file.
 *
 * \param s Handle to finalize.  May be \c NULL.
 * \return ESIO_SUCCESS \c (0) on success or
 *         one of ::esio_error_status on failure.
 */
int esio_finalize(esio_state s);
/*\@}*/


/** \name Opening and closing files */
/*\@{*/

/**
 * Create a new file or overwrite an existing one.
 *
 * \param s Handle to use.
 * \param file Name of the file to open.
 * \param overwrite If zero, fail if an existing file is detected.
 *                  If nonzero, clobber any existing file.
 *
 * \return ESIO_SUCCESS \c (0) on success or
 *         one of ::esio_error_status on failure.
 */
int esio_file_create(esio_state s, const char *file, int overwrite);

/**
 * Open an existing file.
 *
 * \param s Handle to use.
 * \param file Name of the file to open.
 * \param readwrite If zero, open the file in read-only mode.
 *                  If nonzero, open the file in read-write mode.
 *
 * \return ESIO_SUCCESS \c (0) on success or
 *         one of ::esio_error_status on failure.
 */
int esio_file_open(esio_state s, const char *file, int readwrite);

/**
 * Flush buffers associated with any currently open file.
 *
 * \param s Handle to use.
 *
 * \return ESIO_SUCCESS \c (0) on success or
 *         one of ::esio_error_status on failure.
 */
int esio_file_flush(esio_state s);

/**
 * Close any currently open file.
 * Closing a file automatically flushes all unwritten data.
 *
 * \param s Handle to use.
 *
 * \return ESIO_SUCCESS \c (0) on success or
 *         one of ::esio_error_status on failure.
 */
int esio_file_close(esio_state s);
/*\@}*/


/** \name Manipulating string-valued attributes */
/*\@{*/

/**
 * Set the attribute \c name to be \c value.
 * Any existing attribute will be overwritten.
 *
 * @param s Handle to use.
 * @param name Null-terminated attribute name.
 * @param value Null-terminated attribute value.
 *
 * \return ESIO_SUCCESS \c (0) on success or
 *         one of ::esio_error_status on failure.
 */
int esio_string_set(const esio_state s,
                    const char *name,
                    const char *value);

/**
 * Get the null-terminated string associated with attribute \c name.  The
 * routine allocates sufficient storage to hold the string and returns it.  The
 * caller must <tt>free</tt> the memory to avoid resource leaks.
 *
 * @param s Handle to use.
 * @param name Null-terminated attribute name.
 *
 * @return A newly allocated buffer containing the value on success.
 *         \c NULL on failure.
 */
char* esio_string_get(const esio_state s,
                      const char *name);
/*\@}*/


/** @cond INTERNAL */
#define ESIO_ATTRIBUTE_GEN(TYPE)                     \
int esio_attribute_write_##TYPE(const esio_state s,  \
                                const char *name,    \
                                const TYPE *value);  \
int esio_attribute_read_##TYPE(const esio_state s,   \
                               const char *name,     \
                               TYPE *value);

#define ESIO_ATTRIBUTE_GENV(TYPE)                    \
int esio_attribute_writev_##TYPE(const esio_state s, \
                                 const char *name,   \
                                 const TYPE *value,  \
                                 int ncomponents);   \
int esio_attribute_readv_##TYPE(const esio_state s,  \
                                const char *name,    \
                                TYPE *value,         \
                                int ncomponents);
/** @endcond */

/** \name Manipulating scalar-valued attributes */
/*\@{*/
ESIO_ATTRIBUTE_GEN(double)
ESIO_ATTRIBUTE_GEN(float)
ESIO_ATTRIBUTE_GEN(int)
/*\@}*/

/** \name Manipulating vector-valued attributes */
/*\@{*/
ESIO_ATTRIBUTE_GENV(double)
ESIO_ATTRIBUTE_GENV(float)
ESIO_ATTRIBUTE_GENV(int)
int esio_attribute_sizev(const esio_state s,
                         const char *name,
                         int *ncomponents);
/*\@}*/

#undef ESIO_ATTRIBUTE_GEN
#undef ESIO_ATTRIBUTE_GENV


/** @cond INTERNAL */
#define ESIO_LINE_GEN(TYPE)                                                   \
int esio_line_write_##TYPE(const esio_state s,                                \
                           const char *name,                                  \
                           const TYPE *line,                                  \
                           int aglobal, int astart, int alocal, int astride); \
int esio_line_read_##TYPE(const esio_state s,                                 \
                          const char *name,                                   \
                          TYPE *line,                                         \
                          int aglobal, int astart, int alocal, int astride);

#define ESIO_LINE_GENV(TYPE)                                                  \
int esio_line_writev_##TYPE(const esio_state s,                               \
                            const char *name,                                 \
                            const TYPE *line,                                 \
                            int aglobal, int astart, int alocal, int astride, \
                            int ncomponents);                                 \
int esio_line_readv_##TYPE(const esio_state s,                                \
                           const char *name,                                  \
                           TYPE *line,                                        \
                           int aglobal, int astart, int alocal, int astride,  \
                           int ncomponents);
/** @endcond */

/** \name Manipulating distributed, one-dimensional, scalar-valued data */
/*\@{*/
ESIO_LINE_GEN(double)
ESIO_LINE_GEN(float)
ESIO_LINE_GEN(int)
int esio_line_size(const esio_state s,
                   const char *name,
                   int *aglobal);
/*\@}*/

/** \name Manipulating distributed, one-dimensional, vector-valued data */
/*\@{*/
ESIO_LINE_GENV(double)
ESIO_LINE_GENV(float)
ESIO_LINE_GENV(int)
int esio_line_sizev(const esio_state s,
                    const char *name,
                    int *aglobal,
                    int *ncomponents);
/*\@}*/

#undef ESIO_LINE_GEN
#undef ESIO_LINE_GENV


/** @cond INTERNAL */
#define ESIO_PLANE_GEN(TYPE)                                                  \
int esio_plane_write_##TYPE(const esio_state s,                               \
                            const char *name,                                 \
                            const TYPE *plane,                                \
                            int bglobal, int bstart, int blocal, int bstride, \
                            int aglobal, int astart, int alocal, int astride);\
int esio_plane_read_##TYPE(const esio_state s,                                \
                           const char *name,                                  \
                           TYPE *plane,                                       \
                           int bglobal, int bstart, int blocal, int bstride,  \
                           int aglobal, int astart, int alocal, int astride);

#define ESIO_PLANE_GENV(TYPE)                                                 \
int esio_plane_writev_##TYPE(const esio_state s,                              \
                             const char *name,                                \
                             const TYPE *plane,                               \
                             int bglobal, int bstart, int blocal, int bstride,\
                             int aglobal, int astart, int alocal, int astride,\
                             int ncomponents);                                \
int esio_plane_readv_##TYPE(const esio_state s,                               \
                            const char *name,                                 \
                            TYPE *plane,                                      \
                            int bglobal, int bstart, int blocal, int bstride, \
                            int aglobal, int astart, int alocal, int astride, \
                            int ncomponents);
/** @endcond */

/** \name Manipulating distributed, two-dimensional, scalar-valued data */
/*\@{*/
ESIO_PLANE_GEN(double)
ESIO_PLANE_GEN(float)
ESIO_PLANE_GEN(int)
int esio_plane_size(const esio_state s,
                    const char *name,
                    int *bglobal, int *aglobal);
/*\@}*/

/** \name Manipulating distributed, two-dimensional, vector-valued data */
ESIO_PLANE_GENV(double)
ESIO_PLANE_GENV(float)
ESIO_PLANE_GENV(int)
int esio_plane_sizev(const esio_state s,
                     const char *name,
                     int *bglobal, int *aglobal,
                     int *ncomponents);
/*\@}*/

#undef ESIO_PLANE_GENV
#undef ESIO_PLANE_GEN


/** @cond INTERNAL */
#define ESIO_FIELD_GEN(TYPE)                                                  \
int esio_field_write_##TYPE(const esio_state s,                               \
                            const char *name,                                 \
                            const TYPE *field,                                \
                            int cglobal, int cstart, int clocal, int cstride, \
                            int bglobal, int bstart, int blocal, int bstride, \
                            int aglobal, int astart, int alocal, int astride);\
int esio_field_read_##TYPE(const esio_state s,                                \
                           const char *name,                                  \
                           TYPE *field,                                       \
                           int cglobal, int cstart, int clocal, int cstride,  \
                           int bglobal, int bstart, int blocal, int bstride,  \
                           int aglobal, int astart, int alocal, int astride);

#define ESIO_FIELD_GENV(TYPE)                                                 \
int esio_field_writev_##TYPE(const esio_state s,                              \
                             const char *name,                                \
                             const TYPE *field,                               \
                             int cglobal, int cstart, int clocal, int cstride,\
                             int bglobal, int bstart, int blocal, int bstride,\
                             int aglobal, int astart, int alocal, int astride,\
                             int ncomponents);                                \
int esio_field_readv_##TYPE(const esio_state s,                               \
                            const char *name,                                 \
                            TYPE *field,                                      \
                            int cglobal, int cstart, int clocal, int cstride, \
                            int bglobal, int bstart, int blocal, int bstride, \
                            int aglobal, int astart, int alocal, int astride, \
                            int ncomponents);
/** @endcond */

/** \name Manipulating distributed, three-dimensional, scalar-valued data */
/*\@{*/
ESIO_FIELD_GEN(double)
ESIO_FIELD_GEN(float)
ESIO_FIELD_GEN(int)
int esio_field_size(const esio_state s,
                    const char *name,
                    int *cglobal, int *bglobal, int *aglobal);
/*\@}*/


/** \name Manipulating distributed, three-dimensional, vector-valued data */
/*\@{*/
ESIO_FIELD_GENV(double)
ESIO_FIELD_GENV(float)
ESIO_FIELD_GENV(int)
int esio_field_sizev(const esio_state s,
                     const char *name,
                     int *cglobal, int *bglobal, int *aglobal,
                     int *ncomponents);
/*\@}*/

#undef ESIO_FIELD_GEN
#undef ESIO_FIELD_GENV


/** \name Querying and controlling field layout */
/*\@{*/
int esio_layout_count();
int esio_layout_get(const esio_state s);
int esio_layout_set(esio_state s, int layout_tag);
/*\@}*/

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __ESIO_ESIO_H */
