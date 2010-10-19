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
 * esio.h: Header describing ESIO's public C API
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

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
 * \defgroup Size Querying data extents
 */
/*\@{*/
int esio_point_size(const esio_state s,
                    const char* name);
int esio_point_sizev(const esio_state s,
                     const char* name,
                     int *ncomponents);
int esio_line_size(const esio_state s,
                   const char* name,
                   int *aglobal);
int esio_line_sizev(const esio_state s,
                    const char* name,
                    int *aglobal,
                    int *ncomponents);
int esio_plane_size(const esio_state s,
                    const char* name,
                    int *bglobal, int *aglobal);
int esio_plane_sizev(const esio_state s,
                     const char* name,
                     int *bglobal, int *aglobal,
                     int *ncomponents);
int esio_field_size(const esio_state s,
                    const char* name,
                    int *cglobal, int *bglobal, int *aglobal);
int esio_field_sizev(const esio_state s,
                     const char* name,
                     int *cglobal, int *bglobal, int *aglobal,
                     int *ncomponents);
/*\@}*/

/**
 * \defgroup Point Manipulating real-valued points
 */
/*\@{*/
#define ESIO_POINT_GEN(TYPE)                     \
int esio_point_write_##TYPE(const esio_state s,  \
                            const char* name,    \
                            const TYPE *point);  \
int esio_point_writev_##TYPE(const esio_state s, \
                             const char* name,   \
                             const TYPE *point,  \
                             int ncomponents);   \
int esio_point_read_##TYPE(const esio_state s,   \
                           const char* name,     \
                           TYPE *point);         \
int esio_point_readv_##TYPE(const esio_state s,  \
                            const char* name,    \
                            TYPE *point,         \
                            int ncomponents);
ESIO_POINT_GEN(double)
ESIO_POINT_GEN(float)
ESIO_POINT_GEN(int)
#undef ESIO_POINT_GEN
/*\@}*/

/**
 * \defgroup Line Manipulating real-valued lines
 */
/*\@{*/
#define ESIO_LINE_GEN(TYPE)                                                   \
int esio_line_write_##TYPE(const esio_state s,                                \
                           const char* name,                                  \
                           const TYPE *line,                                  \
                           int aglobal, int astart, int alocal, int astride); \
int esio_line_writev_##TYPE(const esio_state s,                               \
                            const char* name,                                 \
                            const TYPE *line,                                 \
                            int aglobal, int astart, int alocal, int astride, \
                            int ncomponents);                                 \
int esio_line_read_##TYPE(const esio_state s,                                 \
                          const char* name,                                   \
                          TYPE *line,                                         \
                          int aglobal, int astart, int alocal, int astride);  \
int esio_line_readv_##TYPE(const esio_state s,                                \
                           const char* name,                                  \
                           TYPE *line,                                        \
                           int aglobal, int astart, int alocal, int astride,  \
                           int ncomponents);
ESIO_LINE_GEN(double)
ESIO_LINE_GEN(float)
ESIO_LINE_GEN(int)
#undef ESIO_LINE_GEN
/*\@}*/

/**
 * \defgroup Plane Manipulating real-valued planes
 */
/*\@{*/
#define ESIO_PLANE_GEN(TYPE)                                                  \
int esio_plane_write_##TYPE(const esio_state s,                               \
                            const char* name,                                 \
                            const TYPE *plane,                                \
                            int bglobal, int bstart, int blocal, int bstride, \
                            int aglobal, int astart, int alocal, int astride);\
int esio_plane_writev_##TYPE(const esio_state s,                              \
                             const char* name,                                \
                             const TYPE *plane,                               \
                             int bglobal, int bstart, int blocal, int bstride,\
                             int aglobal, int astart, int alocal, int astride,\
                             int ncomponents);                                \
int esio_plane_read_##TYPE(const esio_state s,                                \
                           const char* name,                                  \
                           TYPE *plane,                                       \
                           int bglobal, int bstart, int blocal, int bstride,  \
                           int aglobal, int astart, int alocal, int astride); \
int esio_plane_readv_##TYPE(const esio_state s,                               \
                            const char* name,                                 \
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
 * \defgroup Field Manipulating real-valued fields
 */
/*\@{*/
#define ESIO_FIELD_GEN(TYPE)                                                  \
int esio_field_write_##TYPE(const esio_state s,                               \
                            const char* name,                                 \
                            const TYPE *field,                                \
                            int cglobal, int cstart, int clocal, int cstride, \
                            int bglobal, int bstart, int blocal, int bstride, \
                            int aglobal, int astart, int alocal, int astride);\
int esio_field_writev_##TYPE(const esio_state s,                              \
                             const char* name,                                \
                             const TYPE *field,                               \
                             int cglobal, int cstart, int clocal, int cstride,\
                             int bglobal, int bstart, int blocal, int bstride,\
                             int aglobal, int astart, int alocal, int astride,\
                             int ncomponents);                                \
int esio_field_read_##TYPE(const esio_state s,                                \
                           const char* name,                                  \
                           TYPE *field,                                       \
                           int cglobal, int cstart, int clocal, int cstride,  \
                           int bglobal, int bstart, int blocal, int bstride,  \
                           int aglobal, int astart, int alocal, int astride); \
int esio_field_readv_##TYPE(const esio_state s,                               \
                            const char* name,                                 \
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
