//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// ESIO 0.1.9: ExaScale IO library for turbulence simulation restart files
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "layout.h"

#include <assert.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "error.h"

// *******************************************************************
// PLANE LINE POINT PLANE LINE POINT PLANE LINE POINT PLANE LINE POINT
// *******************************************************************

#define METHODNAME esio_plane_writer
#define OPFUNC     H5Dwrite
#define QUALIFIER  const
#include "x-plane.c"
#undef METHODNAME
#undef OPFUNC
#undef QUALIFIER

#define METHODNAME esio_plane_reader
#define OPFUNC     H5Dread
#define QUALIFIER  /* mutable */
#include "x-plane.c"
#undef METHODNAME
#undef OPFUNC
#undef QUALIFIER

#define METHODNAME esio_line_writer
#define OPFUNC     H5Dwrite
#define QUALIFIER  const
#include "x-line.c"
#undef METHODNAME
#undef OPFUNC
#undef QUALIFIER

#define METHODNAME esio_line_reader
#define OPFUNC     H5Dread
#define QUALIFIER  /* mutable */
#include "x-line.c"
#undef METHODNAME
#undef OPFUNC
#undef QUALIFIER

// ***********************************************************************
// LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0
// ***********************************************************************

hid_t esio_field_layout0_filespace_creator(int cglobal,
                                           int bglobal,
                                           int aglobal)
{
    const hsize_t dims[3] = { cglobal, bglobal, aglobal };
    return H5Screate_simple(3, dims, NULL);
}

herr_t esio_field_layout0_dataset_chunker(hid_t dcpl_id,
                                          int cchunk, int bchunk, int achunk)
{
    const hsize_t chunksizes[3] = { cchunk, bchunk, achunk };
    return H5Pset_chunk(dcpl_id, 3, chunksizes);
}

#define METHODNAME esio_field_layout0_field_writer
#define OPFUNC     H5Dwrite
#define QUALIFIER  const
#include "x-layout0.c"
#undef METHODNAME
#undef OPFUNC
#undef QUALIFIER

#define METHODNAME esio_field_layout0_field_reader
#define OPFUNC     H5Dread
#define QUALIFIER  /* mutable */
#include "x-layout0.c"
#undef METHODNAME
#undef OPFUNC
#undef QUALIFIER

// ***********************************************************************
// LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1
// ***********************************************************************

hid_t esio_field_layout1_filespace_creator(int cglobal,
                                           int bglobal,
                                           int aglobal)
{
    const hsize_t dims[3] = { cglobal, bglobal, aglobal };
    return H5Screate_simple(3, dims, NULL);
}

herr_t esio_field_layout1_dataset_chunker(hid_t dcpl_id,
                                          int cchunk, int bchunk, int achunk)
{
    const hsize_t chunksizes[3] = { cchunk, bchunk, achunk };
    return H5Pset_chunk(dcpl_id, 3, chunksizes);
}

#define METHODNAME esio_field_layout1_field_writer
#define OPFUNC     H5Dwrite
#define QUALIFIER  const
#include "x-layout1.c"
#undef METHODNAME
#undef OPFUNC
#undef QUALIFIER

#define METHODNAME esio_field_layout1_field_reader
#define OPFUNC     H5Dread
#define QUALIFIER  /* mutable */
#include "x-layout1.c"
#undef METHODNAME
#undef OPFUNC
#undef QUALIFIER

// ***********************************************************************
// LAYOUT 2 LAYOUT 2 LAYOUT 2 LAYOUT 2 LAYOUT 2 LAYOUT 2 LAYOUT 2 LAYOUT 2
// ***********************************************************************

hid_t esio_field_layout2_filespace_creator(int cglobal,
                                           int bglobal,
                                           int aglobal)
{
    const hsize_t dims[2] = { cglobal * bglobal, aglobal };
    return H5Screate_simple(2, dims, NULL);
}

herr_t esio_field_layout2_dataset_chunker(hid_t dcpl_id,
                                          int cchunk, int bchunk, int achunk)
{
    const hsize_t chunksizes[2] = { cchunk * bchunk, achunk };
    return H5Pset_chunk(dcpl_id, 2, chunksizes);
}

#define METHODNAME esio_field_layout2_field_writer
#define OPFUNC     H5Dwrite
#define QUALIFIER  const
#include "x-layout2.c"
#undef METHODNAME
#undef OPFUNC
#undef QUALIFIER

#define METHODNAME esio_field_layout2_field_reader
#define OPFUNC     H5Dread
#define QUALIFIER  /* mutable */
#include "x-layout2.c"
#undef METHODNAME
#undef OPFUNC
#undef QUALIFIER
