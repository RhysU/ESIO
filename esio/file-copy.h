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

#ifndef ESIO_FILE_COPY_H
#define ESIO_FILE_COPY_H

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

#endif /* ESIO_FILE_COPY_H */
