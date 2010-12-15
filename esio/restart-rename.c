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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "error.h"
#include "restart-rename.h"

int restart_rename(const char *src_filename,
                   const char *dst_template,
                   int keep_howmany)
{
    if (src_filename == NULL) ESIO_ERROR("src_filename == NULL", ESIO_EFAULT);
    if (dst_template == NULL) ESIO_ERROR("dst_template == NULL", ESIO_EFAULT);
    if (keep_howmany < 0)     ESIO_ERROR("keep_howmany < 0",     ESIO_EINVAL);

    // Ensure we can stat src_filename, which should (mostly) isolate
    // rename(2) ENOENT errors to be related to the destination file.
    struct stat statbuf;
    if (stat(src_filename, &statbuf) < 0) {
        ESIO_ERROR("Error to stat(2)-ing src_filename during restart_rename",
                   ESIO_EFAILED);
    }

    // FIXME Implement
}
