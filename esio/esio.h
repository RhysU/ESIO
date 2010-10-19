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
int esio_field_size(const esio_state s,
                    const char* name,
                    int *cglobal, int *bglobal, int *aglobal);
int esio_vfield_size(const esio_state s,
                     const char* name,
                     int *cglobal, int *bglobal, int *aglobal,
                     int *ncomponents);
/*\@}*/

/**
 * \defgroup Field Manipulating real-valued fields 
 */
/*\@{*/
int esio_field_write_double(const esio_state s,
                            const char* name,
                            const double *field,
                            int cglobal, int cstart, int clocal, int cstride,
                            int bglobal, int bstart, int blocal, int bstride,
                            int aglobal, int astart, int alocal, int astride);
int esio_field_write_float(const esio_state s,
                           const char* name,
                           const float *field,
                           int cglobal, int cstart, int clocal, int cstride,
                           int bglobal, int bstart, int blocal, int bstride,
                           int aglobal, int astart, int alocal, int astride);
int esio_field_read_double(const esio_state s,
                           const char* name,
                           double *field,
                           int cglobal, int cstart, int clocal, int cstride,
                           int bglobal, int bstart, int blocal, int bstride,
                           int aglobal, int astart, int alocal, int astride);
int esio_field_read_float(const esio_state s,
                          const char* name,
                          float *field,
                          int cglobal, int cstart, int clocal, int cstride,
                          int bglobal, int bstart, int blocal, int bstride,
                          int aglobal, int astart, int alocal, int astride);
/*\@}*/

/**
 * \defgroup VField Manipulating real-valued vector fields 
 */
/*\@{*/
int esio_vfield_write_double(const esio_state s,
                             const char* name,
                             const double *field,
                             int cglobal, int cstart, int clocal, int cstride,
                             int bglobal, int bstart, int blocal, int bstride,
                             int aglobal, int astart, int alocal, int astride,
                             int ncomponents);
int esio_vfield_write_float(const esio_state s,
                            const char* name,
                            const float *field,
                            int cglobal, int cstart, int clocal, int cstride,
                            int bglobal, int bstart, int blocal, int bstride,
                            int aglobal, int astart, int alocal, int astride,
                            int ncomponents);
int esio_vfield_read_double(const esio_state s,
                            const char* name,
                            double *field,
                            int cglobal, int cstart, int clocal, int cstride,
                            int bglobal, int bstart, int blocal, int bstride,
                            int aglobal, int astart, int alocal, int astride,
                            int ncomponents);
int esio_vfield_read_float(const esio_state s,
                           const char* name,
                           float *field,
                           int cglobal, int cstart, int clocal, int cstride,
                           int bglobal, int bstart, int blocal, int bstride,
                           int aglobal, int astart, int alocal, int astride,
                           int ncomponents);
/*\@}*/

#endif /* __ESIO_ESIO_H */
