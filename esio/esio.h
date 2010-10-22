//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// esio 0.0.1: ExaScale IO library for turbulence simulation restart files
// 
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

/* TODO Document C API */

#ifndef __ESIO_ESIO_H
#define __ESIO_ESIO_H

#include <mpi.h>

typedef struct esio_state_s *esio_state;

/**
 * \defgroup Init Setup and teardown
 */
/*\@{*/
esio_state esio_init(MPI_Comm comm);
int esio_finalize(esio_state s);
/*\@}*/

/**
 * \defgroup Layout Field layout control
 */
/*\@{*/
int esio_layout_count();
int esio_layout_get(const esio_state s);
int esio_layout_set(esio_state s, int layout_tag);

/*\@}*/

/**
 * \defgroup File Opening and closing files
 */
/*\@{*/
int esio_file_create(esio_state s, const char *file, int overwrite);
int esio_file_open(esio_state s, const char *file, int readwrite);
int esio_file_close(esio_state s);
/*\@}*/

/**
 * \defgroup Attribute Manipulating simple attributes
 */
/*\@{*/
int esio_attribute_sizev(const esio_state s,
                         const char *name,
                         int *ncomponents);

#define ESIO_ATTRIBUTE_GEN(TYPE)                     \
int esio_attribute_write_##TYPE(const esio_state s,  \
                                const char *name,    \
                                const TYPE *value);  \
int esio_attribute_writev_##TYPE(const esio_state s, \
                                 const char *name,   \
                                 const TYPE *value,  \
                                 int ncomponents);   \
int esio_attribute_read_##TYPE(const esio_state s,   \
                               const char *name,     \
                               TYPE *value);         \
int esio_attribute_readv_##TYPE(const esio_state s,  \
                                const char *name,    \
                                TYPE *value,         \
                                int ncomponents);
ESIO_ATTRIBUTE_GEN(double)
ESIO_ATTRIBUTE_GEN(float)
ESIO_ATTRIBUTE_GEN(int)
#undef ESIO_ATTRIBUTE_GEN
/*\@}*/

/**
 * \defgroup String Manipulating strings
 */
/*\@{*/
int esio_string_set(const esio_state s,
                    const char *name,
                    const char *value);
char* esio_string_get(const esio_state s,
                      const char *name);
/*\@}*/

/**
 * \defgroup Size Querying distributed data extents
 */
/*\@{*/
int esio_line_size(const esio_state s,
                   const char *name,
                   int *aglobal);
int esio_line_sizev(const esio_state s,
                    const char *name,
                    int *aglobal,
                    int *ncomponents);
int esio_plane_size(const esio_state s,
                    const char *name,
                    int *bglobal, int *aglobal);
int esio_plane_sizev(const esio_state s,
                     const char *name,
                     int *bglobal, int *aglobal,
                     int *ncomponents);
int esio_field_size(const esio_state s,
                    const char *name,
                    int *cglobal, int *bglobal, int *aglobal);
int esio_field_sizev(const esio_state s,
                     const char *name,
                     int *cglobal, int *bglobal, int *aglobal,
                     int *ncomponents);
/*\@}*/


/**
 * \defgroup Line Manipulating distributed 1D lines
 */
/*\@{*/
#define ESIO_LINE_GEN(TYPE)                                                   \
int esio_line_write_##TYPE(const esio_state s,                                \
                           const char *name,                                  \
                           const TYPE *line,                                  \
                           int aglobal, int astart, int alocal, int astride); \
int esio_line_writev_##TYPE(const esio_state s,                               \
                            const char *name,                                 \
                            const TYPE *line,                                 \
                            int aglobal, int astart, int alocal, int astride, \
                            int ncomponents);                                 \
int esio_line_read_##TYPE(const esio_state s,                                 \
                          const char *name,                                   \
                          TYPE *line,                                         \
                          int aglobal, int astart, int alocal, int astride);  \
int esio_line_readv_##TYPE(const esio_state s,                                \
                           const char *name,                                  \
                           TYPE *line,                                        \
                           int aglobal, int astart, int alocal, int astride,  \
                           int ncomponents);
ESIO_LINE_GEN(double)
ESIO_LINE_GEN(float)
ESIO_LINE_GEN(int)
#undef ESIO_LINE_GEN
/*\@}*/

/**
 * \defgroup Plane Manipulating distributed 2D planes
 */
/*\@{*/
#define ESIO_PLANE_GEN(TYPE)                                                  \
int esio_plane_write_##TYPE(const esio_state s,                               \
                            const char *name,                                 \
                            const TYPE *plane,                                \
                            int bglobal, int bstart, int blocal, int bstride, \
                            int aglobal, int astart, int alocal, int astride);\
int esio_plane_writev_##TYPE(const esio_state s,                              \
                             const char *name,                                \
                             const TYPE *plane,                               \
                             int bglobal, int bstart, int blocal, int bstride,\
                             int aglobal, int astart, int alocal, int astride,\
                             int ncomponents);                                \
int esio_plane_read_##TYPE(const esio_state s,                                \
                           const char *name,                                  \
                           TYPE *plane,                                       \
                           int bglobal, int bstart, int blocal, int bstride,  \
                           int aglobal, int astart, int alocal, int astride); \
int esio_plane_readv_##TYPE(const esio_state s,                               \
                            const char *name,                                 \
                            TYPE *plane,                                      \
                            int bglobal, int bstart, int blocal, int bstride, \
                            int aglobal, int astart, int alocal, int astride, \
                            int ncomponents);
ESIO_PLANE_GEN(double)
ESIO_PLANE_GEN(float)
ESIO_PLANE_GEN(int)
#undef ESIO_PLANE_GEN
/*\@}*/

/**
 * \defgroup Field Manipulating distributed 3D fields
 */
/*\@{*/
#define ESIO_FIELD_GEN(TYPE)                                                  \
int esio_field_write_##TYPE(const esio_state s,                               \
                            const char *name,                                 \
                            const TYPE *field,                                \
                            int cglobal, int cstart, int clocal, int cstride, \
                            int bglobal, int bstart, int blocal, int bstride, \
                            int aglobal, int astart, int alocal, int astride);\
int esio_field_writev_##TYPE(const esio_state s,                              \
                             const char *name,                                \
                             const TYPE *field,                               \
                             int cglobal, int cstart, int clocal, int cstride,\
                             int bglobal, int bstart, int blocal, int bstride,\
                             int aglobal, int astart, int alocal, int astride,\
                             int ncomponents);                                \
int esio_field_read_##TYPE(const esio_state s,                                \
                           const char *name,                                  \
                           TYPE *field,                                       \
                           int cglobal, int cstart, int clocal, int cstride,  \
                           int bglobal, int bstart, int blocal, int bstride,  \
                           int aglobal, int astart, int alocal, int astride); \
int esio_field_readv_##TYPE(const esio_state s,                               \
                            const char *name,                                 \
                            TYPE *field,                                      \
                            int cglobal, int cstart, int clocal, int cstride, \
                            int bglobal, int bstart, int blocal, int bstride, \
                            int aglobal, int astart, int alocal, int astride, \
                            int ncomponents);
ESIO_FIELD_GEN(double)
ESIO_FIELD_GEN(float)
ESIO_FIELD_GEN(int)
#undef ESIO_FIELD_GEN
/*\@}*/


#endif /* __ESIO_ESIO_H */
