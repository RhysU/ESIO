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
hid_t esio_field_layout ## NUM ## _filespace_creator(      \
        int cglobal, int bglobal, int aglobal);            \
                                                           \
int esio_field_layout ## NUM ## _field_writer(             \
        hid_t dset_id, const void *field,                  \
        int cglobal, int cstart, int clocal, int cstride,  \
        int bglobal, int bstart, int blocal, int bstride,  \
        int aglobal, int astart, int alocal, int astride,  \
        hid_t type_id);                                    \
                                                           \
int esio_field_layout ## NUM ##_field_reader(              \
        hid_t dset_id, void *field,                        \
        int cglobal, int cstart, int clocal, int cstride,  \
        int bglobal, int bstart, int blocal, int bstride,  \
        int aglobal, int astart, int alocal, int astride,  \
        hid_t type_id);

ESIO_LAYOUT_DECLARATIONS(0)
ESIO_LAYOUT_DECLARATIONS(1)

int esio_plane_writer(
        hid_t dset_id, const void *plane,
        int bglobal, int bstart, int blocal, int bstride,
        int aglobal, int astart, int alocal, int astride,
        hid_t type_id);

int esio_plane_reader(
        hid_t dset_id, void *plane,
        int bglobal, int bstart, int blocal, int bstride,
        int aglobal, int astart, int alocal, int astride,
        hid_t type_id);

int esio_line_writer(
        hid_t dset_id, const void *line,
        int aglobal, int astart, int alocal, int astride,
        hid_t type_id);

int esio_line_reader(
        hid_t dset_id, void *line,
        int aglobal, int astart, int alocal, int astride,
        hid_t type_id);

int esio_point_writer(
        hid_t dset_id, const void *point,
        hid_t type_id);

int esio_point_reader(
        hid_t dset_id, void *point,
        hid_t type_id);

#endif /* __ESIO_LAYOUT_H */
