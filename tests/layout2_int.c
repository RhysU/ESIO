//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// ESIO 0.2.0: ExaScale IO library for turbulence simulation restart files
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

#define REAL              int
#define REAL_H5T          H5T_NATIVE_INT
#define AFFIX(name)       name ## _int
#define LAYOUT_TAG        (2)

#include "layout_template.c"
