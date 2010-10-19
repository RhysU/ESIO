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
 * layout.c: Implementation details for various field layouts
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "error.h"
#include "layout.h"

// *******************************************************************
// PLANE LINE POINT PLANE LINE POINT PLANE LINE POINT PLANE LINE POINT
// *******************************************************************

#define METHODNAME esio_plane_writer
#define OPFUNC     H5Dwrite
#define QUALIFIER  const
#include "plane.c"
#undef METHODNAME
#undef OPFUNC
#undef QUALIFIER

#define METHODNAME esio_plane_reader
#define OPFUNC     H5Dread
#define QUALIFIER  /* mutable */
#include "plane.c"
#undef METHODNAME
#undef OPFUNC
#undef QUALIFIER

#define METHODNAME esio_line_writer
#define OPFUNC     H5Dwrite
#define QUALIFIER  const
#include "line.c"
#undef METHODNAME
#undef OPFUNC
#undef QUALIFIER

#define METHODNAME esio_line_reader
#define OPFUNC     H5Dread
#define QUALIFIER  /* mutable */
#include "line.c"
#undef METHODNAME
#undef OPFUNC
#undef QUALIFIER

// FIXME: Enable
// #define METHODNAME esio_point_writer
// #define OPFUNC     H5Dwrite
// #define QUALIFIER  const
// #include "point.c"
// #undef METHODNAME
// #undef OPFUNC
// #undef QUALIFIER
// 
// #define METHODNAME esio_point_reader
// #define OPFUNC     H5Dread
// #define QUALIFIER  /* mutable */
// #include "point.c"
// #undef METHODNAME
// #undef OPFUNC
// #undef QUALIFIER

// ***********************************************************************
// LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0
// ***********************************************************************

hid_t esio_layout0_filespace_creator(int cglobal, int bglobal, int aglobal)
{
    const hsize_t dims[3] = { cglobal, bglobal, aglobal };
    return H5Screate_simple(3, dims, NULL);
}

#define METHODNAME esio_layout0_field_writer
#define OPFUNC     H5Dwrite
#define QUALIFIER  const
#include "layout0.c"
#undef METHODNAME
#undef OPFUNC
#undef QUALIFIER

#define METHODNAME esio_layout0_field_reader
#define OPFUNC     H5Dread
#define QUALIFIER  /* mutable */
#include "layout0.c"
#undef METHODNAME
#undef OPFUNC
#undef QUALIFIER

// ***********************************************************************
// LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1
// ***********************************************************************

hid_t esio_layout1_filespace_creator(int cglobal, int bglobal, int aglobal)
{
    const hsize_t dims[2] = { cglobal * bglobal, aglobal };
    return H5Screate_simple(2, dims, NULL);
}

#define METHODNAME esio_layout1_field_writer
#define OPFUNC     H5Dwrite
#define QUALIFIER  const
#include "layout1.c"
#undef METHODNAME
#undef OPFUNC
#undef QUALIFIER

#define METHODNAME esio_layout1_field_reader
#define OPFUNC     H5Dread
#define QUALIFIER  /* mutable */
#include "layout1.c"
#undef METHODNAME
#undef OPFUNC
#undef QUALIFIER
