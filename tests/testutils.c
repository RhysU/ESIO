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

#include <assert.h>
#include <errno.h>
#include <libgen.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "testutils.h"

char * create_testfiletemplate(const char *dir, const char * const filename)
{
    // Use a default if dir is NULL or empty
    if (!dir || *dir == 0) dir = P_tmpdir;

    // Determine basename from filename
    assert(filename);
    char *       t1 = strdup(filename);
    char * const t2 = basename(t1);

    // Determine required storage
    size_t len = 0;
    len += strlen(dir);
    len += strlen("/test.");
    len += strlen(t2);
    len += strlen(".XXXXXX");
    len += 1;

    // Allocate required storage
    char * const retval = malloc(len);
    assert(retval);

    // Produce a file template suitable for mkstemp(3)
    strcpy(retval, dir);
    strcat(retval, "/test.");
    strcat(retval, t2);
    strcat(retval, ".XXXXXX");

    // Deallocate temporary storage
    free(t1);

    return retval;
}

char * create_testfilename(const char * const testfiletemplate)
{
    // Create a temporary filename using mkstemp(3)
    // Requires calling mkstemp(3) and deleting the file it creates.
    // Still allows for tempnam(3)-like races, but silences linker.
    // Dicey but not horrible for simple test cases

    assert(testfiletemplate);

    char * retval = strdup(testfiletemplate);

    errno = 0;
    if (retval) {
        close(mkstemp(retval));
        unlink(retval);
        if (errno) {
            free(retval);
            retval = NULL;
        }
    }

    return retval;
}
