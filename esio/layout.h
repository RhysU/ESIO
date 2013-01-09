//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// ESIO 0.1.8: ExaScale IO library for turbulence simulation restart files
// http://red.ices.utexas.edu/projects/esio/
//
// Copyright (C) 2010, 2011, 2012, 2013 The PECOS Development Team
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

#ifndef ESIO_LAYOUT_H
#define ESIO_LAYOUT_H

#include <hdf5.h>

#ifdef __cplusplus
extern "C" {
#endif

//*********************************************************************
// INTERNAL TYPES INTERNAL TYPES INTERNAL TYPES INTERNAL TYPES INTERNAL
//*********************************************************************

typedef hid_t  (*esio_filespace_creator_t)(int, int, int);

typedef herr_t (*esio_dataset_chunker_t)  (hid_t, int, int, int);

typedef int    (*esio_field_writer_t)     (hid_t, hid_t, const void *,
                                           int, int, int, int,
                                           int, int, int, int,
                                           int, int, int, int,
                                           hid_t);

typedef int    (*esio_field_reader_t)     (hid_t, hid_t, void *,
                                           int, int, int, int,
                                           int, int, int, int,
                                           int, int, int, int,
                                           hid_t);

//******************************************************************
// INTERNAL DECLARATIONS INTERNAL DECLARATIONS INTERNAL DECLARATIONS
//******************************************************************

#define ESIO_LAYOUT_DECLARATIONS(NUM)                      \
hid_t esio_field_layout ## NUM ## _filespace_creator(      \
        int cglobal, int bglobal, int aglobal);            \
                                                           \
herr_t esio_field_layout ## NUM ## _dataset_chunker(       \
        hid_t dcpl_id,                                     \
        int cchunk, int bchunk, int achunk);               \
                                                           \
int esio_field_layout ## NUM ## _field_writer(             \
        hid_t plist_id, hid_t dset_id, const void *field,  \
        int cglobal, int cstart, int clocal, int cstride,  \
        int bglobal, int bstart, int blocal, int bstride,  \
        int aglobal, int astart, int alocal, int astride,  \
        hid_t type_id);                                    \
                                                           \
int esio_field_layout ## NUM ##_field_reader(              \
        hid_t plist_id, hid_t dset_id, void *field,        \
        int cglobal, int cstart, int clocal, int cstride,  \
        int bglobal, int bstart, int blocal, int bstride,  \
        int aglobal, int astart, int alocal, int astride,  \
        hid_t type_id);

ESIO_LAYOUT_DECLARATIONS(0)
ESIO_LAYOUT_DECLARATIONS(1)
ESIO_LAYOUT_DECLARATIONS(2)

int esio_plane_writer(
        hid_t plist_id, hid_t dset_id, const void *plane,
        int bglobal, int bstart, int blocal, int bstride,
        int aglobal, int astart, int alocal, int astride,
        hid_t type_id);

int esio_plane_reader(
        hid_t plist_id, hid_t dset_id, void *plane,
        int bglobal, int bstart, int blocal, int bstride,
        int aglobal, int astart, int alocal, int astride,
        hid_t type_id);

int esio_line_writer(
        hid_t plist_id, hid_t dset_id, const void *line,
        int aglobal, int astart, int alocal, int astride,
        hid_t type_id);

int esio_line_reader(
        hid_t plist_id, hid_t dset_id, void *line,
        int aglobal, int astart, int alocal, int astride,
        hid_t type_id);

#ifdef __cplusplus
}
#endif

#endif /* ESIO_LAYOUT_H */
