//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// esio 0.1.4: ExaScale IO library for turbulence simulation restart files
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <errno.h>
#include <fcntl.h>
#include <stddef.h>
#include <stdlib.h>
#include <unistd.h>

#include "binary-io.h"
#include "full-write.h"
#include "safe-read.h"

/* ESIO's error handling */
#include "error.h"
#include "file-copy.h"

int
file_copy(const char *src_filename,
          const char *dest_filename,
          int overwrite,
          int blockuntilsync)
{
    /* Allocate working buffer */
    const int OPSIZE = 32 * 1024;
    char *buf = malloc(OPSIZE);
    if (!buf) {
        ESIO_ERROR("Unable to allocate buffer for copy operation",
                   ESIO_ENOMEM);
    }

    /* Open up source and destination files */
    const int src = open(src_filename, O_RDONLY | O_BINARY);
    if (src < 0) {
        ESIO_ERROR("Error opening copy source for reading", ESIO_EFAILED);
    }

    /* Open a new file or truncate an existing one based on overwrite */
    int dst;
    if (overwrite) {
        dst = open(dest_filename,
                   O_BINARY | O_CREAT | O_TRUNC | O_WRONLY, 0600);
        if (dst < 0) {
            ESIO_ERROR("Error overwriting copy destination file for writing",
                        ESIO_EFAILED);
        }
    } else {
        dst = open(dest_filename,
                   O_BINARY | O_CREAT | O_EXCL | O_WRONLY, 0600);
        if (dst < 0) {
            ESIO_ERROR("Error opening new copy destination file for writing",
                       ESIO_EFAILED);
        }
    }

    /* Copy file contents in chunks of size OPSIZE */
    while (1) {
        const size_t n = safe_read(src, buf, OPSIZE);
        if (n == SAFE_READ_ERROR) {
            ESIO_ERROR("Error reading source file", ESIO_EFAILED);
        }
        if (n == 0) {
            break;
        }

        if (full_write(dst, buf, n) < n) {
            ESIO_ERROR("Error writing destination file", ESIO_EFAILED);
        }
    }

    /* Done with temporary buffer */
    free(buf);

    /* If requested, ensure the destination file has hit the device */
    if (blockuntilsync) {
        if (fsync(dst) < 0) {
            ESIO_ERROR("Error sync-ing destination file", ESIO_EFAILED);
        }
    }

    /* Close source and destination file */
    if (close(dst) < 0) {
        ESIO_ERROR("Error closing destination file", ESIO_EFAILED);
    }
    if (close(src) < 0) {
        ESIO_ERROR("Error closing source file", ESIO_EFAILED);
    }

    return ESIO_SUCCESS;
}
