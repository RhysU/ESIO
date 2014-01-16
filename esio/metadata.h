//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// ESIO 0.1.9: ExaScale IO library for turbulence simulation restart files
// http://red.ices.utexas.edu/projects/esio/
//
// Copyright (C) 2010-2014 The PECOS Development Team
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

#ifndef ESIO_METADATA_H
#define ESIO_METADATA_H

#include <hdf5.h>

#ifdef __cplusplus
extern "C" {
#endif

//*********************************************************************
// INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL PROTOTYPES INTERNAL
//*********************************************************************

int esio_type_ncomponents(hid_t type_id);

hid_t esio_type_arrayify(hid_t type_id, int ncomponents);

int esio_field_metadata_write(hid_t loc_id, const char *name,
                              int layout_index,
                              int cglobal, int bglobal, int aglobal,
                              int ncomponents);

int esio_field_metadata_read(hid_t loc_id, const char *name,
                             int *layout_index,
                             int *cglobal, int *bglobal, int *aglobal,
                             int *ncomponents);

int esio_plane_metadata_write(hid_t loc_id, const char *name,
                              int bglobal, int aglobal,
                              int ncomponents);

int esio_plane_metadata_read(hid_t loc_id, const char *name,
                             int *bglobal, int *aglobal,
                             int *ncomponents);

int esio_line_metadata_write(hid_t loc_id, const char *name,
                             int aglobal,
                             int ncomponents);

int esio_line_metadata_read(hid_t loc_id, const char *name,
                            int *aglobal,
                            int *ncomponents);

#ifdef __cplusplus
}
#endif

#endif /* ESIO_METADATA_H */
