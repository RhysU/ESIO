/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PSDNS Development Team
 *
 * Please see https://wiki.ices.utexas.edu/PSDNS for more information.
 *
 * This file is part of the esio library.
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
 * layout0_float.c: single precision unit tests for ESIO
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#define TEST_REAL              float
#define TEST_ESIO_FIELD_WRITE  esio_field_write_float
#define TEST_ESIO_FIELD_READ   esio_field_read_float
#define TEST_ESIO_VFIELD_WRITE esio_vfield_write_float
#define TEST_ESIO_VFIELD_READ  esio_vfield_read_float
#define TEST_H5T               H5T_NATIVE_FLOAT
#define LAYOUT_TAG             (0)

#include "layout_template.c"
