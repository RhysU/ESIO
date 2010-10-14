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

esio_state esio_init(MPI_Comm comm);
int esio_finalize(esio_state s);

int esio_file_create(esio_state s, const char *file, int overwrite);
int esio_file_open(esio_state s, const char *file, int readwrite);
int esio_file_close(esio_state s);

int esio_field_size(esio_state s,
                    const char* name,
                    int *nc, int *nb, int *na);

int esio_field_write_double(esio_state s,
                            const char* name,
                            const double *field,
                            int nc, int cst, int csz,
                            int nb, int bst, int bsz,
                            int na, int ast, int asz);
int esio_field_write_float(esio_state s,
                           const char* name,
                           const float *field,
                           int nc, int cst, int csz,
                           int nb, int bst, int bsz,
                           int na, int ast, int asz);

int esio_field_read_double(esio_state s,
                           const char* name,
                           double *field,
                           int nc, int cst, int csz,
                           int nb, int bst, int bsz,
                           int na, int ast, int asz);
int esio_field_read_float(esio_state s,
                          const char* name,
                          float *field,
                          int nc, int cst, int csz,
                          int nb, int bst, int bsz,
                          int na, int ast, int asz);


#endif /* __ESIO_ESIO_H */
