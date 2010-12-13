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

#ifndef __ESIO_VISIBILITY_H__
#define __ESIO_VISIBILITY_H__

// Content and ideas taken from http://gcc.gnu.org/wiki/Visibility

// Generic helper definitions for shared library support
#if defined _WIN32 || defined __CYGWIN__
  #define ESIO_HELPER_DLL_IMPORT __declspec(dllimport)
  #define ESIO_HELPER_DLL_EXPORT __declspec(dllexport)
  #define ESIO_HELPER_DLL_LOCAL
#else
  #if __GNUC__ >= 4
    #define ESIO_HELPER_DLL_IMPORT __attribute__ ((visibility("default")))
    #define ESIO_HELPER_DLL_EXPORT __attribute__ ((visibility("default")))
    #define ESIO_HELPER_DLL_LOCAL  __attribute__ ((visibility("hidden")))
  #else
    #define ESIO_HELPER_DLL_IMPORT
    #define ESIO_HELPER_DLL_EXPORT
    #define ESIO_HELPER_DLL_LOCAL
  #endif
#endif

// Now we use the generic helper definitions above to define ESIO_API and
// ESIO_LOCAL.  ESIO_API is used for the public API symbols.  It either DLL
// imports or DLL exports (or does nothing for static build).  ESIO_LOCAL is
// used for non-API symbols.

#ifdef ESIO_DLL // defined if ESIO is compiled as a DLL
  #ifdef ESIO_DLL_EXPORTS // defined if building the ESIO DLL versus using it
    #define ESIO_API ESIO_HELPER_DLL_EXPORT
  #else
    #define ESIO_API ESIO_HELPER_DLL_IMPORT
  #endif // ESIO_DLL_EXPORTS
  #define ESIO_LOCAL ESIO_HELPER_DLL_LOCAL
#else // ESIO_DLL is not defined: this means ESIO is a static lib.
  #define ESIO_API
  #define ESIO_LOCAL
#endif // ESIO_DLL

#endif /* __ESIO_VISIBILITY_H__ */
