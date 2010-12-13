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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* Specification.  */
#include "file-copy.h"

#include <errno.h>
#include <fcntl.h>
#include <stddef.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>

#include "safe-read.h"
#include "full-write.h"
#include "binary-io.h"

/* ESIO, not gnulib, error handling */
#include "error.h"

enum { IO_SIZE = 32 * 1024 };

int
file_copy (const char *src_filename,
           const char *dest_filename,
           int overwrite)
{
    int src_fd;
    struct stat statbuf;
    int mode;
    int dest_fd;
    char *buf = malloc (IO_SIZE);
    if (!buf) {
        ESIO_ERROR("Unable to allocate buffer for copy operation",
                   ESIO_ENOMEM);
    }

    src_fd = open (src_filename, O_RDONLY | O_BINARY);
    if (src_fd < 0 || fstat (src_fd, &statbuf) < 0) {
        ESIO_ERROR("Error opening copy source file for reading",
                   ESIO_EFAILED);
    }

    mode = statbuf.st_mode & 07777;

    /* Open a new file or truncate an existing one based on overwrite */
    if (overwrite) {
        dest_fd = open (dest_filename,
                        O_WRONLY | O_CREAT | O_TRUNC | O_BINARY, 0600);
        if (dest_fd < 0) {
            ESIO_ERROR("Error overwriting copy destination file for writing",
                        ESIO_EFAILED);
        }
    } else {
        dest_fd = open (dest_filename,
                        O_WRONLY | O_CREAT | O_EXCL | O_BINARY, 0600);
        if (dest_fd < 0) {
            ESIO_ERROR("Error opening new copy destination file for writing",
                       ESIO_EFAILED);
        }
    }

    /* Copy the file contents.  */
    for (;;) {
        size_t n_read = safe_read (src_fd, buf, IO_SIZE);
        if (n_read == SAFE_READ_ERROR) {
            ESIO_ERROR("Error reading source file", ESIO_EFAILED);
        }
        if (n_read == 0) {
            break;
        }

        if (full_write (dest_fd, buf, n_read) < n_read) {
            ESIO_ERROR("Error writing destination file", ESIO_EFAILED);
        }
    }

    free (buf);

    if (close (dest_fd) < 0) {
        ESIO_ERROR("Error closing destination file", ESIO_EFAILED);
    }
    if (close (src_fd) < 0) {
        ESIO_ERROR("Error closing source file", ESIO_EFAILED);
    }

#if HAVE_CHOWN
    /* Preserve the owner and group.  */
    chown (dest_filename, statbuf.st_uid, statbuf.st_gid);
#endif

    /* Preserve the access permissions.  */
    chmod (dest_filename, mode);

    return ESIO_SUCCESS;
}
