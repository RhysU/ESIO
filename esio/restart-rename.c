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
#include <errno.h>
#include <libgen.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "error.h"
#include "restart-rename.h"

static inline int max(int a, int b)
{
    return a > b ? a : b;
}

static inline int min(int a, int b)
{
    return a < b ? a : b;
}

int
snprintf_realloc(char **str, size_t *size, const char *format, ...)
{
    int retval, needed;
    va_list ap;
    va_start(ap, format);
    while (   0     <= (retval = vsnprintf(*str, *size, format, ap)) // Error?
           && *size <  (needed = retval + 1)) {                      // Space?
        va_end(ap);
        *size *= 2;                                                  // Space?
        if (*size < needed) *size = needed;                          // Space!
        char *p = realloc(*str, *size);                              // Alloc
        if (p) {
            *str = p;
        } else {
            free(*str);
            *str  = NULL;
            *size = 0;
            return -1;
        }
        va_start(ap, format);
    }
    va_end(ap);
    return retval;
}

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

// Thread local storage of filter_tmpl, if possible
#ifndef TLS
# define TLS
#endif
static TLS const char *filter_tmpl = NULL;

// Used as an argument to scandir within restart_rename
static int filter(const struct dirent *d)
{
    return restart_nextindex(filter_tmpl, d->d_name, 0);
}

// Used as an argument to scandir within restart_rename
static int compar(const struct dirent **a,
                  const struct dirent **b)
{
    return strverscmp((*a)->d_name, (*b)->d_name);
}

int restart_rename(const char *src_filepath,
                   const char *dst_template,
                   int keep_howmany)
{
    char errmsg[256] = ""; // Used to provide fairly extensive error messages

    if (src_filepath == NULL) ESIO_ERROR("src_filepath == NULL", ESIO_EFAULT);
    if (dst_template == NULL) ESIO_ERROR("dst_template == NULL", ESIO_EFAULT);
    if (keep_howmany < 1)     ESIO_ERROR("keep_howmany < 1",     ESIO_EINVAL);

    // Ensure we can stat src_filepath, which should (mostly) isolate
    // rename(2) ENOENT errors to be related to the destination file.
    struct stat statbuf;
    if (stat(src_filepath, &statbuf) < 0) {
        snprintf(errmsg, sizeof(errmsg),
                 "Error stat(2)-ing src_filepath '%s' during restart_rename",
                 src_filepath);
        ESIO_ERROR(errmsg, ESIO_EFAILED);
    }

    // Split dst_template into dst_dirname/dst_basename
    // Need some auxiliary memory since dirname/basename mutate their argument
    const size_t tmpl_len = strlen(dst_template);
    char *buffer = malloc(/* NB */ 3 * (tmpl_len + 1));
    if (!buffer) ESIO_ERROR("Unable to allocate buffer", ESIO_ENOMEM);
    strcpy(buffer,                dst_template);
    strcpy(buffer + tmpl_len + 1, dst_template);
    const char * tmpl_dirname  = dirname(buffer);
    const char * tmpl_basename = basename(buffer + tmpl_len + 1);

    // Obtain prefix/suffix pointers from splitting tmpl_basename at ('#')*
    // Additionally, count the number of '#' characters within tmpl_basename
    char *prefix = buffer + 2*(tmpl_len + 1);
    strcpy(prefix, tmpl_basename);
    char *suffix  = strchr(prefix, '#');
    if (!suffix) {
        snprintf(errmsg, sizeof(errmsg),
                 "dst_template '%s' must contain at least one '#'",
                 dst_template);
        free(buffer);
        ESIO_ERROR(errmsg, ESIO_EINVAL);
    }
    int ndigits = 0;
    while (*suffix == '#') {
        ++ndigits;
        *suffix++ = 0;
    }

    // Ensure tmpl_basename had only a single sequence of '#' characters
    if (strchr(suffix, '#')) {
        snprintf(errmsg, sizeof(errmsg),
                 "dst_template '%s' cannot contain multiple nonadjacent '#'s",
                 dst_template);
        free(buffer);
        ESIO_ERROR(errmsg, ESIO_EINVAL);
    }

    // Ensure ndigits is not bigger than INT_MAX would make reasonable
    ndigits = min(ndigits, (int) ceil(log(INT_MAX - 1)/log(10)));

    // Ensure ndigits is big enough to comfortably handle keep_howmany choice
    ndigits = max(ndigits, (keep_howmany == 1)
                           ? 1
                           : (int) ceil(log(keep_howmany-1)/log(10)));

    // Scan the directory tmpl_dirname looking for things matching tmpl_basename
    // Sort is according to strverscmp which does what we need here
    filter_tmpl = tmpl_basename; // Thread local, if possible
    struct dirent **namelist;
    int n = scandir(tmpl_dirname, &namelist, filter, compar);
    if (n < 0) {
        snprintf(errmsg, sizeof(errmsg),
                 "Error invoking scandir on '%s' within restart_rename",
                 tmpl_dirname);
        free(buffer);
        ESIO_ERROR(errmsg, ESIO_EFAILED);
    }

    // Buffers in which we'll build source and destination paths
    // These will be intelligently realloc-ed by snprintf_realloc.
    char *srcbuf = NULL, *dstbuf = NULL;
    size_t srclen = 0, dstlen = 0;

    // Process found entries in reverse order.
    // Process means renaming them to be the next index value where appropriate.
    while (n--) {

        // Compute the target index value for the current entry.
        const int next
            = restart_nextindex(tmpl_basename, namelist[n]->d_name, -1);

        // Skip processing any entry which is malformed or out-of-range
        if (next <= 0 || next >= keep_howmany) {
            free(namelist[n]);
            continue;
        }

        // Construct the source pathname for any in-range entry
        if (0 > snprintf_realloc(&srcbuf, &srclen, "%s/%s",
                                 tmpl_dirname, namelist[n]->d_name)) {
            snprintf(errmsg, sizeof(errmsg), "Unable to form srcbuf for '%s'",
                     namelist[n]->d_name);
            if (srcbuf) free(srcbuf);
            if (dstbuf) free(dstbuf);
            do { free(namelist[n]); } while (n--);
            free(namelist);
            free(buffer);
            ESIO_ERROR(errmsg, ESIO_ENOMEM);
        }

        // Construct the destination pathname from the template details
        if (0 > snprintf_realloc(&dstbuf, &dstlen, "%s/%s%0*d%s",
                                 tmpl_dirname, prefix, ndigits, next, suffix)) {
            snprintf(errmsg, sizeof(errmsg), "Unable to form dstbuf for '%s'",
                     namelist[n]->d_name);
            if (srcbuf) free(srcbuf);
            if (dstbuf) free(dstbuf);
            do { free(namelist[n]); } while (n--);
            free(namelist);
            free(buffer);
            ESIO_ERROR(errmsg, ESIO_ENOMEM);
        }

        // Rename the entry from srcbuf to dstbuf using constructed names
        if (rename(srcbuf, dstbuf)) {
            snprintf(errmsg, sizeof(errmsg),
                     "Error renaming '%s' to '%s'", srcbuf, dstbuf);
            if (srcbuf) free(srcbuf);
            if (dstbuf) free(dstbuf);
            do { free(namelist[n]); } while (n--);
            free(namelist);
            free(buffer);
            ESIO_ERROR(errmsg, ESIO_EFAILED);
        }

        // Free current namelist entry
        free(namelist[n]);
    }
    free(namelist);
    if (srcbuf) free(srcbuf);

    // Build the destination for the src_filepath rename
    if (0 > snprintf_realloc(&dstbuf, &dstlen, "%s/%s%0*d%s",
                             tmpl_dirname, prefix, ndigits, 0, suffix)) {
        snprintf(errmsg, sizeof(errmsg), "Unable to form dstbuf for '%s'",
                 src_filepath);
        if (dstbuf) free(dstbuf);
        free(buffer);
        ESIO_ERROR(errmsg, ESIO_ENOMEM);
    }

    // Finally, perform the rename of src_filepath to match dst_template
    if (rename(src_filepath, dstbuf)) {
        snprintf(errmsg, sizeof(errmsg),
                 "Error renaming '%s' to '%s'", src_filepath, dstbuf);
        if (dstbuf) free(dstbuf);
        free(buffer);
        ESIO_ERROR(errmsg, ESIO_EFAILED);
    }

    // Clean up
    free(dstbuf);
    free(buffer);

    return ESIO_SUCCESS;
}
