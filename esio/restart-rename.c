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

#include <assert.h>
#include <ctype.h>
#include <dirent.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "error.h"
#include "restart-rename.h"

int restart_nextindex(const char *tmpl,
                      const char *name,
                      const int errval)
{
    assert(tmpl);
    assert(name);

    // Advance both until the first hash is encountered
    int i = 0;
    while (tmpl[i] && name[i] && tmpl[i] == name[i]) ++i;
    if (tmpl[i] == 0)      return errval; // Usage error
    if (tmpl[i] != '#')    return  0;     // Mismatch
    if (!isdigit(name[i])) return  0;     // Mismatch and/or leading sign

    // Advance template until its end, looking for the final hash
    int j = i, k = i + 1;
    while (tmpl[k]) {
        if (tmpl[k] == '#') j = k;
        ++k;
    }

    // Advance name until its end
    int l = i + 1;
    while (name[l]) ++l;

    // Scan both backwards until the final hash is encountered
    while (k > j && l > i && tmpl[k] == name[l]) {
        --k;
        --l;
    }
    if (tmpl[k] != '#') return 0; // Mismatch

    // Attempt to read a decimal unsigned long from name[i, l]
    char *endptr = NULL;
    const unsigned long curr = strtoul(name + i, &endptr, 10);
    if (endptr != name + l + 1) return 0;      // Mismatch
    if (curr > INT_MAX - 1)     return errval; // Overflow

    // Sanity check that template contained only a single hash sequence
    while (i != j) {
        if (tmpl[i] != '#') return errval; // Usage error
        ++i;
    }

    // Increment current and return.  Overflow prevented above.
    return ((int) curr) + 1;
}

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

    // Sanity check arguments
    // Ensure src_filename present
    // Split dst_template into dst_dirname/dst_basename
    // Ensure dst_dirname exists
    // Ensure dst_basename well-formed (only 1 sequence of # characters)
    //    Calculate offsets for pre, ###, post
    // Calculate ndigits base on maximum of
    //    1) ceil(log(keep_howmany)/log(10.0))
    //    2) length of # character sequence
    // Convert dst_basename to dst_basename_printf by replacing ### with %0NNd
    //     where NN = ndigits

    // Use scandir to obtain a sorted list of current files matching pattern:
    //    1) Use strverscmp to sort versions per versionsort manpage suggestion
    //    2) Write a custom filter that:
    //        2a) Checks for a name match
    //        2b) Checks for a regular file
}
