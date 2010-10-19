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
 * metadata.h: (INTERNAL) Implementation details around ESIO's metadata
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __ESIO_METADATA_H
#define __ESIO_METADATA_H

//*********************************************************************
// INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL
//*********************************************************************

#include <hdf5.h>

int esio_type_ncomponents(hid_t type_id);

hid_t esio_type_arrayify(hid_t type_id, int ncomponents);

herr_t esio_field_metadata_write(hid_t loc_id, const char *name,
                                 int layout_tag,
                                 int cglobal, int bglobal, int aglobal,
                                 int ncomponents);

herr_t esio_field_metadata_read(hid_t loc_id, const char *name,
                                int *layout_tag,
                                int *cglobal, int *bglobal, int *aglobal,
                                int *ncomponents);

herr_t esio_plane_metadata_write(hid_t loc_id, const char *name,
                                 int bglobal, int aglobal,
                                 int ncomponents);

herr_t esio_plane_metadata_read(hid_t loc_id, const char *name,
                                int *bglobal, int *aglobal,
                                int *ncomponents);

herr_t esio_line_metadata_write(hid_t loc_id, const char *name,
                                int aglobal,
                                int ncomponents);

herr_t esio_line_metadata_read(hid_t loc_id, const char *name,
                               int *aglobal,
                               int *ncomponents);

herr_t esio_point_metadata_write(hid_t loc_id, const char *name,
                                 int ncomponents);

herr_t esio_point_metadata_read(hid_t loc_id, const char *name,
                                int *ncomponents);

#endif /* __ESIO_METADATA_H */
