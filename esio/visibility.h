//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// esio 0.1.5: ExaScale IO library for turbulence simulation restart files
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

#ifndef __ESIO_VISIBILITY_H__
#define __ESIO_VISIBILITY_H__

/** \cond INTERNAL */

// Content and ideas modified from http://gcc.gnu.org/wiki/Visibility

// Generic helper definitions for shared library support
#if __GNUC__ >= 4
#  define ESIO_HELPER_SHARED_IMPORT __attribute__ ((visibility("default")))
#  define ESIO_HELPER_SHARED_EXPORT __attribute__ ((visibility("default")))
#  define ESIO_HELPER_SHARED_LOCAL  __attribute__ ((visibility("hidden")))
#else
#  define ESIO_HELPER_SHARED_IMPORT
#  define ESIO_HELPER_SHARED_EXPORT
#  define ESIO_HELPER_SHARED_LOCAL
#endif

// Now we use the generic helper definitions above to define ESIO_API and
// ESIO_LOCAL.  ESIO_API is used for the public API symbols.  It either imports
// or exports the relevant symbol.  ESIO_LOCAL is used for non-API symbols.

#ifdef ESIO_SHARED_EXPORTS // defined if building ESIO versus just using it
#  define ESIO_API ESIO_HELPER_SHARED_EXPORT
#else
#  define ESIO_API ESIO_HELPER_SHARED_IMPORT
#endif // ESIO_SHARED_EXPORTS
#define ESIO_LOCAL ESIO_HELPER_SHARED_LOCAL

/** \endcond */

#endif /* __ESIO_VISIBILITY_H__ */
