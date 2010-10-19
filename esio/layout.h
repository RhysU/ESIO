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
 * layout.h: (INTERNAL) Declarations for the various field layout options
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __ESIO_LAYOUT_H
#define __ESIO_LAYOUT_H

#include <hdf5.h>

//*********************************************************************
// INTERNAL TYPES INTERNAL TYPES INTERNAL TYPES INTERNAL TYPES INTERNAL
//*********************************************************************

typedef hid_t (*esio_filespace_creator_t)(int, int, int);

typedef int   (*esio_field_writer_t)     (hid_t, const void *,
                                          int, int, int, int,
                                          int, int, int, int,
                                          int, int, int, int,
                                          hid_t);

typedef int   (*esio_field_reader_t)     (hid_t, void *,
                                          int, int, int, int,
                                          int, int, int, int,
                                          int, int, int, int,
                                          hid_t);

//******************************************************************
// INTERNAL DECLARATIONS INTERNAL DECLARATIONS INTERNAL DECLARATIONS
//******************************************************************

#define ESIO_LAYOUT_DECLARATIONS(NUM)                      \
hid_t esio_layout ## NUM ## _filespace_creator(            \
        int cglobal, int bglobal, int aglobal);            \
                                                           \
int esio_layout ## NUM ## _field_writer(                   \
        hid_t dset_id, const void *field,                  \
        int cglobal, int cstart, int clocal, int cstride,  \
        int bglobal, int bstart, int blocal, int bstride,  \
        int aglobal, int astart, int alocal, int astride,  \
        hid_t type_id);                                    \
                                                           \
int esio_layout ## NUM ##_field_reader(                    \
        hid_t dset_id, void *field,                        \
        int cglobal, int cstart, int clocal, int cstride,  \
        int bglobal, int bstart, int blocal, int bstride,  \
        int aglobal, int astart, int alocal, int astride,  \
        hid_t type_id);

ESIO_LAYOUT_DECLARATIONS(0)
ESIO_LAYOUT_DECLARATIONS(1)

#endif /* __ESIO_LAYOUT_H */
