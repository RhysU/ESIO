//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// ESIO 0.1.8: ExaScale IO library for turbulence simulation restart files
// http://red.ices.utexas.edu/projects/esio/
//
// Copyright (C) 2010, 2011, 2012, 2013 The PECOS Development Team
//
// This file is part of ESIO.
//
// ESIO is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 2.1 of the License, or
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

/******************************************************************************
 * Functionality adopted from the GNU Scientific Library.
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * Relicensed as LGPL by permission of Brian Gough on 19 Feb 2010:
 *   "If anyone wants to use the error macro definitions or
 *   error handler from GSL under the LGPL I am fine with that."
 * See http://www.mail-archive.com/gsl-discuss@sourceware.org/msg00764.html
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "error.h"
#include <stdlib.h>
#include <mpi.h>

esio_error_handler_t * esio_error_handler = NULL;

static void
no_error_handler(const char *reason,
                 const char *file,
                 int line,
                 int esio_errno);

void
esio_error(const char * reason,
           const char * file,
           int line,
           int esio_errno)
{
    if (esio_error_handler) {
        (*esio_error_handler) (reason, file, line, esio_errno);
        return ;
    }

    esio_stream_printf ("ERROR", file, line, reason);

    fflush (stdout);
    fprintf (stderr, "Default esio error handler invoked.\n");
    fflush (stderr);

    MPI_Abort (MPI_COMM_WORLD, 1);
}

esio_error_handler_t *
esio_set_error_handler(esio_error_handler_t * new_handler)
{
    esio_error_handler_t * previous_handler = esio_error_handler;
    esio_error_handler = new_handler;
    return previous_handler;
}


esio_error_handler_t *
esio_set_error_handler_off(void)
{
    esio_error_handler_t * previous_handler = esio_error_handler;
    esio_error_handler = no_error_handler;
    return previous_handler;
}

#ifdef __INTEL_COMPILER
#pragma warning(push,disable:869)
#endif
static void
no_error_handler(const char *reason,
                 const char *file,
                 int line,
                 int esio_errno)
{
    (void) reason;     /* unused */
    (void) file;       /* unused */
    (void) line;       /* unused */
    (void) esio_errno; /* unused */
    return;            /* do nothing */
}
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

FILE * esio_stream = NULL ;
esio_stream_handler_t * esio_stream_handler = NULL;

void
esio_stream_printf(const char *label,
                   const char *file,
                   int line,
                   const char *reason)
{
    if (esio_stream == NULL) {
        esio_stream = stderr;
    }
    if (esio_stream_handler) {
        (*esio_stream_handler) (label, file, line, reason);
        return;
    }
    fprintf(esio_stream,
            "esio: %s:%d: %s: %s\n", file, line, label, reason);

}

esio_stream_handler_t *
esio_set_stream_handler(esio_stream_handler_t * new_handler)
{
    esio_stream_handler_t * previous_handler = esio_stream_handler;
    esio_stream_handler = new_handler;
    return previous_handler;
}

FILE *
esio_set_stream(FILE * new_stream)
{
    FILE * previous_stream;
    if (esio_stream == NULL) {
        esio_stream = stderr;
    }
    previous_stream = esio_stream;
    esio_stream = new_stream;
    return previous_stream;
}

const char *
esio_strerror(const int esio_errno)
{
    switch (esio_errno) {
    case ESIO_SUCCESS:
        return "success" ;
    case ESIO_EFAULT:
        return "invalid pointer" ;
    case ESIO_EINVAL:
        return "invalid argument supplied by user" ;
    case ESIO_EFAILED:
        return "generic failure" ;
    case ESIO_ESANITY:
        return "sanity check failed - shouldn't happen" ;
    case ESIO_ENOMEM:
        return "malloc failed" ;
    case ESIO_NOTFOUND:
        return "object not found" ;
    default:
        return "unknown error code" ;
    }
}
