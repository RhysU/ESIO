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
                                          int, int, int,
                                          int, int, int,
                                          int, int, int,
                                          hid_t, size_t);

//******************************************************************
// INTERNAL DECLARATIONS INTERNAL DECLARATIONS INTERNAL DECLARATIONS
//******************************************************************

hid_t esio_layout0_filespace_creator(int na, int nb, int nc);

int esio_layout0_field_writer(hid_t dset_id, const void *field,
                              int na, int ast, int asz,
                              int nb, int bst, int bsz,
                              int nc, int cst, int csz,
                              hid_t type_id, size_t type_size);

#endif /* __ESIO_LAYOUT_H */
