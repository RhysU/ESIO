/* Copying of files modified from gnulib's copy-file module.
 * Modified to fit into ESIO's error handling framework.
 *
 * Copyright (C) 2001-2003, 2009-2010 Free Software Foundation, Inc.
 * Written by Bruno Haible <haible@clisp.cons.org>, 2001.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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
int file_copy (const char *src_filename,
               const char *dest_filename,
               int overwrite,
               int blockuntilsync);

#ifdef __cplusplus
}
#endif

#endif /* __ESIO_FILE_COPY_H */
