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

//****************************************************************
// DATA STRUCTURES DATA STRUCTURES DATA STRUCTURES DATA STRUCTURES
//****************************************************************

struct field_details {
    MPI_Comm comm;
    esio_handle h;
    int cglobal, cstart, clocal;
    int bglobal, bstart, blocal;
    int aglobal, astart, alocal;
    int dims[3], coords[3], rank;
    long bytes;
    void *data;
};

struct plane_details {
    MPI_Comm comm;
    esio_handle h;
    int bglobal, bstart, blocal;
    int aglobal, astart, alocal;
    int dims[2], coords[2], rank;
    long bytes;
    void *data;
};

struct line_details {
    MPI_Comm comm;
    esio_handle h;
    int aglobal, astart, alocal;
    int dims[1], coords[1], rank;
    long bytes;
    void *data;
};

struct details {
    int world_rank;
    int world_size;
    int verbose;
    int repeat;
    int ncomponents;
    int nfields, nplanes, nlines;
    struct field_details *f;
    struct plane_details *p;
    struct line_details  *l;
};

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

static int field_initialize(const struct details *d, struct field_details *f);

static int plane_initialize(const struct details *d, struct plane_details *p);

static int line_initialize( const struct details *d, struct line_details  *l);

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
    {"verbose",     'v', 0,       0, "produce verbose output", -1 },
    {"repeat",      'r', "count", 0, "number of repetitions", -1 },
    {"ncomponents", 'n', "count", 0, "number of components",   0 },
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

// Parse a single option following Argp semantics
static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
    // Get the input argument from argp_parse.
    struct details *d = state->input;

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
                    case 0:  n = &d->nfields; break;
                    case 1:  n = &d->nplanes; break;
                    case 2:  n = &d->nlines;  break;
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
            if (d->nfields <= 0 && d->nfields <= 0 && d->nplanes <= 0) {
                argp_error(state,
                           "At least one of %s must be strictly positive",
                           args_doc);
            }
            break;

        case 'v':
            d->verbose = 1;
            break;

        case 'r':
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%d %c", &d->repeat, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "repeat option is malformed: '%s'", arg);
            }
            if (d->repeat < 1) {
                argp_failure(state, EX_USAGE, 0,
                        "repeat value %d must be strictly positive",
                        d->repeat);
            }
            break;

        case 'n':
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%d %c",
                            &d->ncomponents, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "ncomponents option is malformed: '%s'", arg);
            }
            if (d->ncomponents < 1) {
                argp_failure(state, EX_USAGE, 0,
                        "ncomponents value %d must be strictly positive",
                        d->ncomponents);
            }
            break;

        case 'f':
            if (d->f->cglobal || d->f->bglobal || d->f->aglobal ) {
                argp_error(state, "only one of --field-{memory,global}"
                           " may be specified");
            }
            if (arg) {
                d->f->bytes = from_human_readable_byte_count(arg);
                if (d->f->bytes < 1)
                    argp_failure(state, EX_USAGE, 0,
                            "field-memory option is malformed: '%s'", arg);
            }
            break;

        case 'p':
            if (   d->p->bglobal || d->p->aglobal) {
                argp_error(state, "only one of --plane-{memory,global}"
                           " may be specified");
            }
            if (arg) {
                d->p->bytes = from_human_readable_byte_count(arg);
                if (d->p->bytes < 1)
                    argp_failure(state, EX_USAGE, 0,
                            "plane-memory option is malformed: '%s'", arg);
            }
            break;

        case 'l':
            if (d->l->aglobal) {
                argp_error(state, "only one of --line-{memory,global}"
                           " may be specified");
            }
            if (arg) {
                d->l->bytes = from_human_readable_byte_count(arg);
                if (d->l->bytes < 1)
                    argp_failure(state, EX_USAGE, 0,
                            "line-memory option is malformed: '%s'", arg);
            }
            break;

        case 'F':
            errno = 0;
            if (3 != sscanf(arg ? arg : "", "%d x %d x %d %c",
                            &d->f->dims[0], &d->f->dims[1],
                            &d->f->dims[2], &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "field-dims option is malformed: '%s'", arg);
            }
            if (d->f->dims[0] < 0 || d->f->dims[1] < 0 || d->f->dims[2] < 0) {
                argp_failure(state, EX_USAGE, 0,
                        "field-dims values %dx%dx%d must be nonnegative",
                        d->f->dims[0], d->f->dims[1], d->f->dims[2]);
            }
            break;


        case 'P':
            errno = 0;
            if (2 != sscanf(arg ? arg : "", "%d x %d %c",
                            &d->p->dims[0], &d->p->dims[1], &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "plane-dims option is malformed: '%s'", arg);
            }
            if (d->p->dims[0] < 0 || d->p->dims[1] < 0) {
                argp_failure(state, EX_USAGE, 0,
                        "plane-dims values %dx%d must be nonnegative",
                        d->p->dims[0], d->p->dims[1]);
            }
            break;


        case 'L':
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%d %c",
                            &d->l->dims[0], &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "line-global option is malformed: '%s'", arg);
            }
            if (d->l->dims[0] < 0) {
                argp_failure(state, EX_USAGE, 0,
                        "line-dims value '%d' must be nonnegative",
                        d->l->dims[0]);
            }
            break;

        case FIELD_GLOBAL:
            if (d->f->bytes) {
                argp_error(state, "only one of --field-{memory,global}"
                           " may be specified");
            }
            errno = 0;
            if (3 != sscanf(arg ? arg : "", "%d x %d x %d %c",
                            &d->f->cglobal, &d->f->bglobal,
                            &d->f->aglobal, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "field-global option not of form CxBxA: '%s'", arg);
            }
            if (d->f->cglobal < 1 || d->f->bglobal < 1 || d->f->aglobal < 1) {
                argp_failure(state, EX_USAGE, 0,
                        "field-global values %dx%dx%d must be strictly positive",
                        d->f->cglobal, d->f->bglobal, d->f->aglobal);
            }
            break;


        case PLANE_GLOBAL:
            if (d->p->bytes) {
                argp_error(state, "only one of --plane-{memory,global}"
                           " may be specified");
            }
            errno = 0;
            if (2 != sscanf(arg ? arg : "", "%d x %d %c",
                            &d->p->bglobal, &d->p->aglobal, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "plane-global option not of form BxA: '%s'", arg);
            }
            if (d->p->bglobal < 1 || d->p->aglobal < 1) {
                argp_failure(state, EX_USAGE, 0,
                        "plane-global values %dx%d must be strictly positive",
                        d->p->bglobal, d->p->aglobal);
            }
            break;

        case LINE_GLOBAL:
            if (d->l->bytes) {
                argp_error(state, "only one of --line-{memory,global}"
                           " may be specified");
            }
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%d %c",
                             &d->l->aglobal, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "line-global option is malformed: '%s'", arg);
            }
            if (d->l->aglobal < 0) {
                argp_failure(state, EX_USAGE, 0,
                        "line-global value %d must be strictly positive",
                        d->l->aglobal);
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

// Rank-dependent output streams established in main().
static FILE *rankout, *rankerr;

int main(int argc, char *argv[])
{
    // Initialize default argument storage and default values
    struct details d;        memset(&d, 0, sizeof(struct details));
    struct field_details f;  memset(&f, 0, sizeof(struct field_details));
    struct plane_details p;  memset(&p, 0, sizeof(struct plane_details));
    struct line_details  l;  memset(&l, 0, sizeof(struct line_details));
    d.repeat = 1;
    d.ncomponents = 1;
    d.f = &f;
    d.p = &p;
    d.l = &l;
    f.comm = p.comm = l.comm = MPI_COMM_NULL;

    // Initialize/finalize MPI
    MPI_Init(&argc, &argv);
    atexit((void (*) ()) MPI_Finalize);
    MPI_Comm_size(MPI_COMM_WORLD, &d.world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &d.world_rank);

    // Establish rank-dependent output streams
    errno   = 0;
    if (d.world_rank == 0) {
        rankout = stdout;
        rankerr = stderr;
    } else {
        rankout = fopen("/dev/null", "w");
        rankerr = rankout;
    }
    if (!rankout) {
        perror("Unable to open rank-dependent output streams"),
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Parse command line arguments using MPI-savvy argp extension
    mpi_argp_parse(d.world_rank, &argp, argc, argv, 0, 0, &d);

    // Initialize the field, plane, and line problems
    if (d.nfields) field_initialize(&d, &f);
    if (d.nplanes) plane_initialize(&d, &p);
    if (d.nlines)  line_initialize( &d, &l);

    // DEBUG: Dump arguments
    printf("verbose:      %d\n", d.verbose);
    printf("ncomponents:  %d\n", d.ncomponents);
    printf("f.global: %d x %d x %d\n", f.cglobal, f.bglobal, f.aglobal);
    printf("f.dims:   %d x %d x %d\n", f.dims[0], f.dims[1], f.dims[2]);
    printf("f.bytes:  %ld\n", f.bytes);
    {
        double coeff;
        const char *units;
        to_human_readable_byte_count(f.bytes, 0, &coeff, &units);
        printf("f.bytes:  %.1f %s\n", coeff, units);

    }
    printf("p.global: %d x %d\n", p.bglobal, p.aglobal);
    printf("p.dims:   %d x %d\n", p.dims[0], p.dims[1]);
    printf("p.bytes:  %ld\n", p.bytes);
    printf("l.global: %d\n", l.aglobal);
    printf("l.dims:   %d\n", l.dims[0]);
    printf("l.bytes:  %ld\n", l.bytes);

/*     esio_handle h = esio_handle_initialize(MPI_COMM_WORLD); */
/*     esio_file_create(h, "data.h5", 1 |+ overwrite +|); */

/*     esio_string_set(h, "program", argv[0]); */

/*     int version = 1; */
/*     esio_attribute_write_int(h, "version", &version); */

/*     double example[2] = { 2.0 * world_rank, 2.0 * world_rank + 1 }; */
/*     esio_l.write_double(h, "example", example, */
/*                            2*world_size, 2*world_rank, 2, 1); */

/*     esio_file_close(h); */
/*     esio_handle_finalize(h); */

    // Free ESIO handles
    esio_handle_finalize(f.h);
    esio_handle_finalize(p.h);
    esio_handle_finalize(l.h);

    // Free MPI communicators
    if (f.comm != MPI_COMM_NULL) MPI_Comm_free(&f.comm);
    if (p.comm != MPI_COMM_NULL) MPI_Comm_free(&p.comm);
    if (l.comm != MPI_COMM_NULL) MPI_Comm_free(&l.comm);

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

static int field_initialize(const struct details *d, struct field_details *f)
{
    // Establish MPI topology for field problem
    ESIO_MPICHKQ(MPI_Dims_create(d->world_size, 3, f->dims));
    int periods[3] = { 0, 0, 0 };
    ESIO_MPICHKQ(MPI_Cart_create(
                MPI_COMM_WORLD, 3, f->dims, periods, 1, &f->comm));
    ESIO_MPICHKQ(MPI_Comm_rank(f->comm, &f->rank));
    ESIO_MPICHKQ(MPI_Cart_coords(f->comm, f->rank, 3, f->coords));

    // Compute global problem size and local portion
    if (f->bytes) {
        const double vectors = d->world_size * f->bytes / d->ncomponents;
        f->cglobal = f->bglobal = f->aglobal = ceil(cbrt(vectors));
    } else {
        f->bytes = f->cglobal * f->bglobal * f->aglobal * d->ncomponents
                 * sizeof(double);
    }
    f->clocal = local(f->cglobal, d->world_size, f->coords[0], 1);
    f->cstart = start(f->cglobal, d->world_size, f->coords[0], 1);
    f->blocal = local(f->bglobal, d->world_size, f->coords[1], 0);
    f->bstart = start(f->bglobal, d->world_size, f->coords[1], 0);
    f->alocal = local(f->aglobal, d->world_size, f->coords[2], 1);
    f->astart = start(f->aglobal, d->world_size, f->coords[2], 1);

    // Initialize ESIO to carry out problem
    f->h = esio_handle_initialize(f->comm);
    esio_field_establish(f->h, f->cglobal, f->clocal, f->cstart,
                               f->bglobal, f->blocal, f->bstart,
                               f->aglobal, f->alocal, f->astart);

    return ESIO_SUCCESS;
}

static int plane_initialize(const struct details *d, struct plane_details *p)
{
    // Establish MPI topology for plane problem
    ESIO_MPICHKQ(MPI_Dims_create(d->world_size, 2, p->dims));
    int periods[2] = { 0, 0 };
    ESIO_MPICHKQ(MPI_Cart_create(
                MPI_COMM_WORLD, 2, p->dims, periods, 1, &p->comm));
    ESIO_MPICHKQ(MPI_Comm_rank(p->comm, &p->rank));
    ESIO_MPICHKQ(MPI_Cart_coords(p->comm, p->rank, 2, p->coords));

    // Compute global problem size and local portion
    if (p->bytes) {
        const double vectors = d->world_size * p->bytes / d->ncomponents;
        p->bglobal = p->aglobal = ceil(sqrt(vectors));
    } else {
        p->bytes = p->bglobal * p->aglobal * d->ncomponents * sizeof(double);
    }
    p->blocal = local(p->bglobal, d->world_size, p->coords[0], 0);
    p->bstart = start(p->bglobal, d->world_size, p->coords[0], 0);
    p->alocal = local(p->aglobal, d->world_size, p->coords[1], 1);
    p->astart = start(p->aglobal, d->world_size, p->coords[1], 1);

    // Initialize ESIO to carry out problem
    p->h = esio_handle_initialize(p->comm);
    esio_plane_establish(p->h, p->bglobal, p->blocal, p->bstart,
                               p->aglobal, p->alocal, p->astart);

    return ESIO_SUCCESS;
}

static int line_initialize(const struct details *d, struct line_details  *l)
{
    // Establish MPI topology for line problem
    ESIO_MPICHKQ(MPI_Dims_create(d->world_size, 1, l->dims));
    int periods[1] = { 0 };
    ESIO_MPICHKQ(MPI_Cart_create(
                MPI_COMM_WORLD, 1, l->dims, periods, 1, &l->comm));
    ESIO_MPICHKQ(MPI_Comm_rank(l->comm, &l->rank));
    ESIO_MPICHKQ(MPI_Cart_coords(l->comm, l->rank, 1, l->coords));

    // Compute global problem size and local portion
    if (l->bytes) {
        const double vectors = d->world_size * l->bytes / d->ncomponents;
        l->aglobal = ceil(vectors);
    } else {
        l->bytes = l->aglobal * d->ncomponents * sizeof(double);
    }
    l->alocal = local(l->aglobal, d->world_size, l->coords[0], 1);
    l->astart = start(l->aglobal, d->world_size, l->coords[0], 1);

    // Initialize ESIO to carry out problem
    l->h = esio_handle_initialize(l->comm);
    esio_line_establish(l->h, l->aglobal, l->alocal, l->astart);

    return ESIO_SUCCESS;
}
