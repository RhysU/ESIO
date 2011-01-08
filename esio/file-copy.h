//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// esio 0.1.1: ExaScale IO library for turbulence simulation restart files
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

#ifndef __ESIO_FILE_COPY_H
#define __ESIO_FILE_COPY_H

//****************************************************************
// INTERNAL INTERNAL INTERNAL INTERNAL INTERNAL INTERNAL INTERNAL
//****************************************************************

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Copy a regular file from \c src_filename to \c dest_filename.
 * Neither permissions nor ownership details are preserved.
 *
 * \param src_filename Source filename
 * \param dest_filename Destination filename
 * \param overwrite If zero, fail if an existing file is detected.
 *                  If nonzero, clobber any existing file.
 * \param blockuntilsync If nonzero, block until the device reports
 *                       that all data has been flushed cleanly.
 *
 * \return Either ESIO_SUCCESS \c (0) or one of ::esio_status on failure.
 */
int file_copy(const char *src_filename,
              const char *dest_filename,
              int overwrite,
              int blockuntilsync);

#ifdef __cplusplus
}
#endif

#endif /* __ESIO_FILE_COPY_H */
