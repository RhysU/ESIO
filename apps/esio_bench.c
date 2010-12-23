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
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sysexits.h>

#include "argp.h"
#include "mpi_argp.h"

#include <mpi.h>
#include <esio/error.h>
#include <esio/esio.h>

//*************************************************************
// STATIC PROTOTYPES STATIC PROTOTYPES STATIC PROTOTYPES STATIC 
//*************************************************************

static void trim(char *a);

static inline int min(int a, int b);

static void to_human_readable_byte_count(long bytes,
                                         int si,
                                         double *coeff,
                                         const char **units);

static long from_human_readable_byte_count(const char *str);

static int local(int nglobal, int nranks, int rank, int extralow);

static int start(int nglobal, int nranks, int rank, int extralow);

//*******************************************************************
// ARGP DETAILS: http://www.gnu.org/s/libc/manual/html_node/Argp.html
//*******************************************************************

const char *argp_program_version      = "esio_bench " PACKAGE_VERSION;
const char *argp_program_bug_address  = PACKAGE_BUGREPORT;
static const char doc[]               =
"Simulate and benchmark ESIO-based application restart write operations."
"\v"
"Write NFIELDS fields, NPLANES planes, and NLINES lines with the "
"specified problem sizes and parallel decompositions using ESIO's "
"restart writing capabilities.  Timing information is collected over "
"one or more iterations.\n";
static const char args_doc[]          = "NFIELDS NPLANES NLINES";

enum {
    FIELD_GLOBAL = 255 /* isascii */,
    PLANE_GLOBAL,
    LINE_GLOBAL
};

static struct argp_option options[] = {
    {"verbose",     'v', 0,       0, "Produce verbose output", -1 },
    {"repeat",      'r', "count", 0, "Number of repetitions", -1 },
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
    int repeat;
    int ncomponents;
    int nfields, nplanes, nlines;
    int field_cglobal, field_bglobal, field_aglobal;
    int plane_bglobal, plane_aglobal;
    int line_aglobal;
    int field_dims[3], plane_dims[2], line_dims[1];
    long field_bytes, plane_bytes, line_bytes;
};

// Parse a single option following Argp semantics
static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
    // Get the input argument from argp_parse.
    struct arguments *arguments = state->input;

    // Trim any leading/trailing whitespace from arg
    if (arg) trim(arg);

    // Want to ensure we consume the entire argument for many options
    // Many sscanf calls provide an extra sentinel %c dumping into &ignore
    char ignore = '\0';

    switch (key) {
        case ARGP_KEY_ARG:
            {
                int *n;
                switch (state->arg_num) {
                    case 0:  n = &arguments->nfields; break;
                    case 1:  n = &arguments->nplanes; break;
                    case 2:  n = &arguments->nlines;  break;
                    default: argp_usage(state);
                }
                errno = 0;
                if (1 != sscanf(arg ? arg : "", "%d %c", n, &ignore)) {
                    argp_failure(state, EX_USAGE, errno,
                            "argument %d is malformed: '%s'",
                            state->arg_num + 1, arg);
                }
                if (*n < 0) {
                    argp_failure(state, EX_USAGE, 0,
                            "argument %d value %d must be nonnegative",
                            state->arg_num + 1, *n);
                }
            }
            break;

        case ARGP_KEY_END:
            if (state->arg_num < 3) {
                argp_usage(state);
            }
            if (    arguments->nfields <= 0
                 && arguments->nfields <= 0
                 && arguments->nplanes <= 0) {
                argp_error(state,
                           "At least one of %s must be strictly positive",
                           args_doc);
            }
            break;

        case 'v':
            arguments->verbose = 1;
            break;

        case 'r':
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%d %c",
                            &arguments->repeat, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "repeat option is malformed: '%s'", arg);
            }
            if (arguments->repeat < 1) {
                argp_failure(state, EX_USAGE, 0,
                        "repeat value %d must be strictly positive",
                        arguments->repeat);
            }
            break;

        case 'n':
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%d %c",
                            &arguments->ncomponents, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "ncomponents option is malformed: '%s'", arg);
            }
            if (arguments->ncomponents < 1) {
                argp_failure(state, EX_USAGE, 0,
                        "ncomponents value %d must be strictly positive",
                        arguments->ncomponents);
            }
            break;

        case 'f':
            if (   arguments->field_cglobal
                || arguments->field_bglobal
                || arguments->field_aglobal ) {
                argp_error(state, "only one of --field-{memory,global}"
                           " may be specified");
            }
            if (arg) {
                arguments->field_bytes = from_human_readable_byte_count(arg);
                if (arguments->field_bytes < 1)
                    argp_failure(state, EX_USAGE, 0,
                            "field-memory option is malformed: '%s'", arg);
            }
            break;

        case 'p':
            if (   arguments->plane_bglobal || arguments->plane_aglobal) {
                argp_error(state, "only one of --plane-{memory,global}"
                           " may be specified");
            }
            if (arg) {
                arguments->plane_bytes = from_human_readable_byte_count(arg);
                if (arguments->plane_bytes < 1)
                    argp_failure(state, EX_USAGE, 0,
                            "plane-memory option is malformed: '%s'", arg);
            }
            break;

        case 'l':
            if (arguments->line_aglobal) {
                argp_error(state, "only one of --line-{memory,global}"
                           " may be specified");
            }
            if (arg) {
                arguments->line_bytes = from_human_readable_byte_count(arg);
                if (arguments->line_bytes < 1)
                    argp_failure(state, EX_USAGE, 0,
                            "line-memory option is malformed: '%s'", arg);
            }
            break;

        case 'F':
            errno = 0;
            if (3 != sscanf(arg ? arg : "", "%d x %d x %d %c",
                            &arguments->field_dims[0],
                            &arguments->field_dims[1],
                            &arguments->field_dims[2], &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "field-dims option is malformed: '%s'", arg);
            }
            if (    arguments->field_dims[0] < 0
                 || arguments->field_dims[1] < 0
                 || arguments->field_dims[2] < 0) {
                argp_failure(state, EX_USAGE, 0,
                        "field-dims values %dx%dx%d must be nonnegative",
                        arguments->field_dims[0],
                        arguments->field_dims[1],
                        arguments->field_dims[2]);
            }
            break;


        case 'P':
            errno = 0;
            if (2 != sscanf(arg ? arg : "", "%d x %d %c",
                            &arguments->plane_dims[0],
                            &arguments->plane_dims[1], &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "plane-dims option is malformed: '%s'", arg);
            }
            if (arguments->plane_dims[0] < 0 || arguments->plane_dims[1] < 0) {
                argp_failure(state, EX_USAGE, 0,
                        "plane-dims values %dx%d must be nonnegative",
                        arguments->plane_dims[0], arguments->plane_dims[1]);
            }
            break;


        case 'L':
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%d %c",
                            &arguments->line_dims[0], &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "line-global option is malformed: '%s'", arg);
            }
            if (arguments->line_dims[0] < 0) {
                argp_failure(state, EX_USAGE, 0,
                        "line-dims value '%d' must be nonnegative",
                        arguments->line_dims[0]);
            }
            break;

        case FIELD_GLOBAL:
            if (arguments->field_bytes) {
                argp_error(state, "only one of --field-{memory,global}"
                           " may be specified");
            }
            errno = 0;
            if (3 != sscanf(arg ? arg : "", "%d x %d x %d %c",
                            &arguments->field_cglobal,
                            &arguments->field_bglobal,
                            &arguments->field_aglobal, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "field-global option not of form CxBxA: '%s'", arg);
            }
            if (    arguments->field_cglobal < 1
                 || arguments->field_bglobal < 1
                 || arguments->field_aglobal < 1) {
                argp_failure(state, EX_USAGE, 0,
                        "field-global values %dx%dx%d must be strictly positive",
                        arguments->field_cglobal,
                        arguments->field_bglobal,
                        arguments->field_aglobal);
            }
            break;


        case PLANE_GLOBAL:
            if (arguments->plane_bytes) {
                argp_error(state, "only one of --plane-{memory,global}"
                           " may be specified");
            }
            errno = 0;
            if (2 != sscanf(arg ? arg : "", "%d x %d %c",
                            &arguments->plane_bglobal,
                            &arguments->plane_aglobal, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "plane-global option not of form BxA: '%s'", arg);
            }
            if (arguments->plane_bglobal < 1 || arguments->plane_aglobal < 1) {
                argp_failure(state, EX_USAGE, 0,
                        "plane-global values %dx%d must be strictly positive",
                        arguments->plane_bglobal, arguments->plane_aglobal);
            }
            break;

        case LINE_GLOBAL:
            if (arguments->line_bytes) {
                argp_error(state, "only one of --line-{memory,global}"
                           " may be specified");
            }
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%d %c",
                            &arguments->line_aglobal, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "line-global option is malformed: '%s'", arg);
            }
            if (arguments->line_aglobal < 0) {
                argp_failure(state, EX_USAGE, 0,
                        "line-global value %d must be strictly positive",
                        arguments->line_aglobal);
            }
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
    arguments.repeat = 1;
    arguments.ncomponents = 1;

    // Parse command line arguments
    mpi_argp_parse(world_rank, &argp, argc, argv, 0, 0, &arguments);

    // Determine MPI_Dims_create-based parallel decompositions
    ESIO_MPICHKQ(MPI_Dims_create(world_size, 3, arguments.field_dims));
    ESIO_MPICHKQ(MPI_Dims_create(world_size, 2, arguments.plane_dims));
    ESIO_MPICHKQ(MPI_Dims_create(world_size, 1, arguments.line_dims));

    // Create MPI communicators for each decomposition and obtain ranks
    MPI_Comm field_comm, plane_comm, line_comm;
    {
        int periods[3] = { 0, 0, 0 };
        ESIO_MPICHKQ(MPI_Cart_create(MPI_COMM_WORLD, 3, arguments.field_dims,
                                     periods, 1, &field_comm));
        ESIO_MPICHKQ(MPI_Cart_create(MPI_COMM_WORLD, 2, arguments.plane_dims,
                                     periods, 1, &plane_comm));
        ESIO_MPICHKQ(MPI_Cart_create(MPI_COMM_WORLD, 1, arguments.line_dims,
                                     periods, 1, &line_comm));
    }
    // Determine this rank's position within each of those decompositions
    int field_rank, plane_rank, line_rank;
    ESIO_MPICHKQ(MPI_Comm_rank(field_comm, &field_rank));
    ESIO_MPICHKQ(MPI_Comm_rank(plane_comm, &plane_rank));
    ESIO_MPICHKQ(MPI_Comm_rank(line_comm,  &line_rank));
    int field_coords[3], plane_coords[2], line_coords[1];
    ESIO_MPICHKQ(MPI_Cart_coords(field_comm, field_rank, 3, field_coords));
    ESIO_MPICHKQ(MPI_Cart_coords(plane_comm, plane_rank, 2, plane_coords));
    ESIO_MPICHKQ(MPI_Cart_coords(line_comm,  line_rank,  1, line_coords));

    // Create ESIO handles for each of those decompositions.
    // Ideally a single communicator would suffice, but the MPI
    // topology details seem to be 1-to-1 with communicators.
    esio_handle field_h = esio_handle_initialize(field_comm);
    esio_handle plane_h = esio_handle_initialize(plane_comm);
    esio_handle line_h  = esio_handle_initialize(line_comm);
    assert(field_h);
    assert(plane_h);
    assert(line_h);

    // Establish field problem details and determine local portion
    if (arguments.field_bytes) {
        double field_vectors = world_size * arguments.field_bytes
                             / arguments.ncomponents;
        arguments.field_cglobal = arguments.field_bglobal // Uniform
                                = arguments.field_aglobal // Uniform
                                = ceil(cbrt(field_vectors));
    } else {
        arguments.field_bytes = arguments.field_cglobal
                              * arguments.field_bglobal
                              * arguments.field_aglobal
                              * arguments.ncomponents
                              * sizeof(double);
    }
    const int field_clocal = local(arguments.field_cglobal, world_size,
                                   field_coords[0], 1);
    const int field_cstart = start(arguments.field_cglobal, world_size,
                                   field_coords[0], 1);
    const int field_blocal = local(arguments.field_bglobal, world_size,
                                   field_coords[1], 0);
    const int field_bstart = start(arguments.field_bglobal, world_size,
                                   field_coords[1], 0);
    const int field_alocal = local(arguments.field_aglobal, world_size,
                                   field_coords[2], 1);
    const int field_astart = start(arguments.field_aglobal, world_size,
                                   field_coords[2], 1);
    esio_field_establish(
            field_h, arguments.field_cglobal, field_clocal, field_cstart,
                     arguments.field_bglobal, field_blocal, field_bstart,
                     arguments.field_aglobal, field_alocal, field_astart);

    // Establish plane problem details and determine local portion
    if (arguments.plane_bytes) {
        double plane_vectors = world_size * arguments.plane_bytes
                             / arguments.ncomponents;
        arguments.plane_bglobal = arguments.plane_aglobal // Uniform
                                = ceil(sqrt(plane_vectors));
    } else {
        arguments.plane_bytes = arguments.plane_bglobal
                              * arguments.plane_aglobal
                              * arguments.ncomponents
                              * sizeof(double);
    }
    const int plane_blocal = local(arguments.plane_bglobal, world_size,
                                   plane_coords[0], 0);
    const int plane_bstart = start(arguments.plane_bglobal, world_size,
                                   plane_coords[0], 0);
    const int plane_alocal = local(arguments.plane_aglobal, world_size,
                                   plane_coords[1], 1);
    const int plane_astart = start(arguments.plane_aglobal, world_size,
                                   plane_coords[1], 1);
    esio_plane_establish(
            plane_h, arguments.plane_bglobal, plane_blocal, plane_bstart,
                     arguments.plane_aglobal, plane_alocal, plane_astart);

    // Establish line problem details and determine local portion
    if (arguments.line_bytes) {
        double line_vectors = world_size * arguments.line_bytes
                            / arguments.ncomponents;
        arguments.line_aglobal = ceil(line_vectors);
    } else {
        arguments.line_bytes = arguments.line_aglobal
                             * arguments.ncomponents
                             * sizeof(double);
    }
    const int line_alocal = local(arguments.line_aglobal, world_size,
                                  line_coords[0], 1);
    const int line_astart = start(arguments.line_aglobal, world_size,
                                  line_coords[0], 1);
    esio_line_establish(
            line_h, arguments.line_aglobal, line_alocal, line_astart);


    // DEBUG: Dump arguments
    printf("verbose:      %d\n", arguments.verbose);
    printf("ncomponents:  %d\n", arguments.ncomponents);
    printf("field_global: %d x %d x %d\n", arguments.field_cglobal,
                                           arguments.field_bglobal,
                                           arguments.field_aglobal);
    printf("field_dims:   %d x %d x %d\n", arguments.field_dims[0],
                                           arguments.field_dims[1],
                                           arguments.field_dims[2]);
    printf("field_bytes:  %ld\n", arguments.field_bytes);
    {
        double coeff;
        const char *units;
        to_human_readable_byte_count(arguments.field_bytes, 0, &coeff, &units);
        printf("field_bytes:  %.1f %s\n", coeff, units);

    }
    printf("plane_global: %d x %d\n", arguments.plane_bglobal,
                                      arguments.plane_aglobal);
    printf("plane_dims:   %d x %d\n", arguments.plane_dims[0],
                                      arguments.plane_dims[1]);
    printf("plane_bytes:  %ld\n", arguments.plane_bytes);
    printf("line_global:  %d\n", arguments.line_aglobal);
    printf("line_dims:    %d\n", arguments.line_dims[0]);
    printf("line_bytes:   %ld\n", arguments.line_bytes);

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

    // Free ESIO handles
    esio_handle_finalize(field_h);
    esio_handle_finalize(plane_h);
    esio_handle_finalize(line_h);

    // Free MPI communicators
    if (field_comm != MPI_COMM_NULL) MPI_Comm_free(&field_comm);
    if (plane_comm != MPI_COMM_NULL) MPI_Comm_free(&plane_comm);
    if (line_comm != MPI_COMM_NULL)  MPI_Comm_free(&line_comm);

    return 0;
}


void trim(char *a)
{
    char *b = a;
    while (isspace(*b))   ++b;
    while (*b)            *a++ = *b++;
    *a = '\0';
    while (isspace(*--a)) *a = '\0';
}


static inline int min(int a, int b)
{
    return a < b ? a : b;
}


// Adapted from http://stackoverflow.com/questions/3758606/
// how-to-convert-byte-size-into-human-readable-format-in-java
static void to_human_readable_byte_count(long bytes,
                                         int si,
                                         double *coeff,
                                         const char **units)
{
    // Static lookup table of byte-based SI units
    static const char *suffix[][2] = { { "B",  "B"   },
                                       { "kB", "KiB" },
                                       { "MB", "MiB" },
                                       { "GB", "GiB" },
                                       { "TB", "TiB" },
                                       { "EB", "EiB" },
                                       { "ZB", "ZiB" },
                                       { "YB", "YiB" } };
    int unit = si ? 1000 : 1024;
    int exp = 0;
    if (bytes > 0) {
        exp = min( (int) (log(bytes) / log(unit)),
                   (int) sizeof(suffix) / sizeof(suffix[0]) - 1);
    }
    *coeff = bytes / pow(unit, exp);
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
    const double coeff = strtod(str, &endptr);
    if (errno) return -1;

    // Skip any intermediate white space
    while (isspace(*endptr)) ++endptr;

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

    // Next character must be one of B/empty/whitespace
    switch (toupper(*endptr)) {
        case 'B':
        case ' ':
        case '\t': ++endptr;  break;

        case '\0': goto done;

        default:   return -1;
    }

    // Skip any remaining white space
    while (isspace(*endptr)) ++endptr;

    // Parse error on anything but a null terminator
    if (*endptr) return -1;

done:
    return exp ? coeff * pow(unit, exp / 3) : coeff;
}

// Compute rank's portion of nglobal elements distributed uniformly across
// nranks ranks.  Any "extra" elements are spread across the lower ranks
// when extralow is set.
static int local(int nglobal, int nranks, int rank, int extralow)
{
    if (extralow) {
        return nglobal / nranks + (rank < nglobal % nranks);
    } else {
        return nglobal / nranks + ((nranks - rank - 1) < nglobal % nranks);
    }
}

// Determine the starting, zero-based offset for a particular rank.
// Not particularly efficient, but effective.
static int start(int nglobal, int nranks, int rank, int extralow)
{
    int offset = -1;
    for (int i = 0; i < rank; ++i) {
        offset += local(nglobal, nranks, rank, extralow);
    }
    return offset;
}
