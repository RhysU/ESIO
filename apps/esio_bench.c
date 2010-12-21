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

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "argp.h"
#include "mpi_argp.h"

#include <mpi.h>
#include <esio/esio.h>

// Argp details: http://www.gnu.org/s/libc/manual/html_node/Argp.html
const char *argp_program_version      = "esio_bench " PACKAGE_VERSION;
const char *argp_program_bug_address  = PACKAGE_BUGREPORT;
static const char doc[]               = "ESIO Benchmarking Tool" // Brief
                                        // "\v"                  // Separator
                                        "";                      // Details
static const char args_doc[]          = "";

enum {
    FIELD_GLOBAL = 255 /* isascii */,
    PLANE_GLOBAL,
    LINE_GLOBAL
};

static struct argp_option options[] = {
    {"verbose",     'v', 0,       0, "Produce verbose output", -1 },
    {"ncomponents", 'n', "count", 0, "Number of components",   0 },
    {0, 0, 0, 0,
     "Controlling field problem size (specify at most one)", 0 },
    {"field-memory", 'f',          "bytes", 0, "per-rank field memory", 0 },
    {"field-global", FIELD_GLOBAL, "CxBxA", 0, "field global extents",  0 },
    {0, 0, 0, 0,
     "Controlling plane problem size (specify at most one)", 0 },
    {"plane-memory", 'p',          "bytes", 0, "per-rank plane memory", 0 },
    {"plane-global", PLANE_GLOBAL, "BxA",   0, "plane global extents",  0 },
    {0, 0, 0, 0,
     "Controlling line problem size (specify at most one)", 0 },
    {"line-memory",  'l',          "bytes", 0, "per-rank line memory", 0 },
    {"line-global",  LINE_GLOBAL,  "A",     0, "line global extents",  0 },
    {0, 0, 0, 0,
     "Controlling parallel decomposition per MPI_Dims_create semantics", 0 },
    {"field-dims", 'F', "NCxNBxNA", 0, "field parallel decomposition",   0 },
    {"plane-dims", 'P', "NBxNA",    0, "plane parallel decomposition",   0 },
    {"line-dims",  'L', "NA",       0, "line parallel decomposition",    0 },
    { 0, 0, 0, 0,  0, 0 }
};

// Used by main to communicate with parse_opt
struct arguments {
    int verbose;
    int field_cglobal, field_bglobal, field_aglobal, field_dims[3];
    int plane_bglobal, plane_aglobal, plane_dims[2];
    int line_aglobal, line_dims[1];
    int ncomponents;
};

// Parse a single option following Argp semantics
static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
    (void) arg; // Possibly unused

    // Get the input argument from argp_parse.
    struct arguments *arguments = state->input;

    switch (key) {
        case 'v':
            arguments->verbose = 1;
            break;

        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

static struct argp argp = {
    options, parse_opt, args_doc, doc,
    0 /*children*/, 0 /*help_filter*/, 0 /*argp_domain*/
};

static inline int min(int a, int b)
{
    return a < b ? a : b;
}

// Adapted from http://stackoverflow.com/questions/3758606/
// how-to-convert-byte-size-into-human-readable-format-in-java
static void to_human_readable_byte_count(long bytes,
                                         int si,
                                         double *scaled,
                                         const char **units)
{
    // Static lookup table of byte-based SI units
    // Looking ahead, 'H' is for the unofficial 'hella' = 10^27.  :P
    static const char *suffix[][2] = { { "B",  "B"   },
                                       { "kB", "KiB" },
                                       { "MB", "MiB" },
                                       { "GB", "GiB" },
                                       { "TB", "TiB" },
                                       { "EB", "EiB" },
                                       { "ZB", "ZiB" },
                                       { "YB", "YiB" },
                                       { "HB", "HiB" } };
    const int unit = si ? 1000 : 1024;
    const int exp  = min( (int) (log(bytes) / log(unit)),
                          (int) sizeof(suffix) / sizeof(suffix[0]) );
    *scaled = bytes / pow(unit, exp);
    *units  = suffix[exp][!!si];
}

// Convert strings like the following into byte counts
//    5MB, 5 MB, 5M, 3.7GB, 123b, 456kBytes
// with some amount of forgiveness baked into the parsing.
static long from_human_readable_byte_count(const char *str)
{
    // Parse leading numeric factor
    char *endptr;
    errno = 0;
    const double scaled = strtod(str, &endptr);
    if (errno) return -1;

    // Skip any intermediate white space
    while (isspace(*endptr++)) /* empty */ ;

    // Read off first character which should be an SI prefix
    int exp  = 0;
    int unit = 1024;
    switch (toupper(*endptr)) {
        case 'B':  exp =  0; break;
        case 'K':  exp =  3; break;
        case 'M':  exp =  6; break;
        case 'G':  exp =  9; break;
        case 'T':  exp = 12; break;
        case 'E':  exp = 15; break;
        case 'Z':  exp = 18; break;
        case 'Y':  exp = 21; break;
        case 'H':  exp = 24; break;

        case ' ':
        case '\t':
        case '\0': exp =  0; goto done;

        default:   return -1;
    }
    ++endptr;

    // If an 'i' or 'I' is present use SI factor-of-1000 units
    if (toupper(*endptr) == 'I') {
        ++endptr;
        unit = 1000;
    }

    // Next character should either be 'B' or empty/whitespace
    switch (toupper(*endptr)) {
        case 'B':
        case ' ':
        case '\t':
        case '\0': goto done;

        default:   return -1;
    }
    ++endptr;

    // Next character must be empty/whitespace
    switch (*endptr) {
        case ' ':
        case '\t':
        case '\0': goto done;

        default:   return -1;
    }

done:
    return scaled * pow(unit, exp);
}


int main(int argc, char *argv[])
{
    // Initialize/finalize MPI
    MPI_Init(&argc, &argv);
    atexit((void (*) ()) MPI_Finalize);
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Initialize default argument values
    struct arguments arguments;
    memset(&arguments, 0, sizeof(struct arguments));
    arguments.ncomponents = 1;

    // Parse arguments
    mpi_argp_parse(world_rank, &argp, argc, argv, 0, 0, &arguments);

/*     esio_handle h = esio_handle_initialize(MPI_COMM_WORLD); */
/*     esio_file_create(h, "data.h5", 1 |+ overwrite +|); */

/*     esio_string_set(h, "program", argv[0]); */

/*     int version = 1; */
/*     esio_attribute_write_int(h, "version", &version); */

/*     double example[2] = { 2.0 * world_rank, 2.0 * world_rank + 1 }; */
/*     esio_line_write_double(h, "example", example, */
/*                            2*world_size, 2*world_rank, 2, 1); */

/*     esio_file_close(h); */
/*     esio_handle_finalize(h); */

    return 0;
}
