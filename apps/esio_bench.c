//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// esio 0.1.2: ExaScale IO library for turbulence simulation restart files
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
#include <unistd.h>

#include "argp.h"
#include "mpi_argp.h"

#include <H5public.h>
#include <mpi.h>
#include <esio/error.h>
#include <esio/esio.h>

#ifdef HAVE_GRVY
#include <grvy.h>
#define GRVY_TIMER_BEGIN(id)   grvy_timer_begin(id)
#define GRVY_TIMER_END(id)     grvy_timer_end(id)
#define GRVY_TIMER_FINALIZE()  grvy_timer_finalize()
#define GRVY_TIMER_INIT(id)    grvy_timer_init(id)
#define GRVY_TIMER_RESET()     grvy_timer_reset()
#define GRVY_TIMER_SUMMARIZE() grvy_timer_summarize()
#else
#define GRVY_TIMER_BEGIN(id)
#define GRVY_TIMER_END(id)
#define GRVY_TIMER_FINALIZE()
#define GRVY_TIMER_INIT(id)
#define GRVY_TIMER_RESET()
#define GRVY_TIMER_SUMMARIZE()
#endif

//****************************************************************
// DATA STRUCTURES DATA STRUCTURES DATA STRUCTURES DATA STRUCTURES
//****************************************************************

struct field_details {
    int cglobal, cstart, clocal;
    int bglobal, bstart, blocal;
    int aglobal, astart, alocal;
    int ncomponents;
    int layout, dims[3], coords[3], rank;
    long bytes;
    void *data;
    long stride;
};

struct plane_details {
    int bglobal, bstart, blocal;
    int aglobal, astart, alocal;
    int dims[2], coords[2], rank;
    int ncomponents;
    long bytes;
    void *data;
    long stride;
};

struct line_details {
    int aglobal, astart, alocal;
    int dims[1], coords[1], rank;
    int ncomponents;
    long bytes;
    void *data;
    long stride;
};

struct details {
    int world_rank;
    int world_size;
    size_t typesize;
    esio_handle h;
    int verbose;
    int repeat;
    int nfields, nplanes, nlines;
    char *uncommitted;
    char *dst_template;
    int   retain;
    struct field_details *f;
    struct plane_details *p;
    struct line_details  *l;
};

//*************************************************************
// STATIC PROTOTYPES STATIC PROTOTYPES STATIC PROTOTYPES STATIC
//*************************************************************

static void print_version(FILE *stream, struct argp_state *state);

static void trim(char *a);

static inline int min(int a, int b);

static inline int max(int a, int b);

static void to_human_readable_byte_count(long bytes,
                                         int si,
                                         double *coeff,
                                         const char **units);

static long from_human_readable_byte_count(const char *str);

static int local(int nglobal, int nranks, int rank, int extralow);

static int start(int nglobal, int nranks, int rank, int extralow);

static int global_minmax(long val, long *min, long *max);

static int field_initialize(struct details *d, struct field_details *f);

static int plane_initialize(struct details *d, struct plane_details *p);

static int line_initialize( struct details *d, struct line_details  *l);

static void* malloc_and_fill(struct details *d, const long bytes);

static int field_finalize(struct details *d, struct field_details *f);

static int plane_finalize(struct details *d, struct plane_details *p);

static int line_finalize( struct details *d, struct line_details  *l);

//*******************************************************************
// ARGP DETAILS: http://www.gnu.org/s/libc/manual/html_node/Argp.html
//*******************************************************************

const char *argp_program_version      = "esio_bench " PACKAGE_VERSION;
void (*argp_program_version_hook)(FILE *stream, struct argp_state *state)
                                      = &print_version;
const char *argp_program_bug_address  = PACKAGE_BUGREPORT;
static const char doc[]               =
"Simulate and benchmark ESIO-based application restart write operations."
"\v"
"Write fields, planes, and lines with the specified problem sizes and "
"parallel decompositions using ESIO's restart writing capabilities. "
"With one argument, a restart file named FILENAME is written.  With two "
"arguments, uncommitted restart files are written to UNCOMMITTED and "
"then renamed to match DESTTEMPLATE.  Timing information is collected "
"over one or more iterations.\n"
"\n"
"Options taking a 'bytes' parameter can be given common byte-related "
"units.  For example --field-memory=5G indicates that approximately "
"5 gigabytes of memory should be used on each rank to store field data. "
"SI units like 'Ki' or 'MiB' are also accepted.\n"
;
static const char args_doc[] = "FILENAME\nUNCOMMITTED DESTTEMPLATE";

enum {
    FIELD_LAYOUT = 255 /* !isascii */,
    FIELD_GLOBAL,
    PLANE_GLOBAL,
    LINE_GLOBAL,
    FIELD_NCOMPONENTS,
    PLANE_NCOMPONENTS,
    LINE_NCOMPONENTS,
    NFIELDS,
    NPLANES,
    NLINES
};

static struct argp_option options[] = {
    {"verbose",     'v', 0,       0, "produce verbose output",            0 },
    {"repeat",      'r', "count", 0, "number of repetitions",             0 },
    {"retain",      'R', "count", 0, "number of restart files to retain", 0 },
    {"nfields",      NFIELDS,      "count", 0, "number of fields",    1},
    {"nplanes",      NPLANES,      "count", 0, "number of planes",    1},
    {"nlines",       NLINES,       "count", 0, "number of lines",     1},
    {"field-layout", FIELD_LAYOUT, "index", 0, "field layout to use", 0},
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
     "Controlling number of components used in each problem type", 0 },
    {"field-ncomponents", FIELD_NCOMPONENTS, "count", 0,
            "number of components per field",   0 },
    {"plane-ncomponents", PLANE_NCOMPONENTS, "count", 0,
            "number of components per plane",   0 },
    {"line-ncomponents",  LINE_NCOMPONENTS, "count", 0,
            "number of components per line",   0 },
    {0, 0, 0, 0,
     "Changing the type of data written", 0 },
    {"single",      's', 0,       0, "write single-precision data", 0 },
    {"double",      'd', 0,       0, "write double-precision data", 0 },
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
            switch (state->arg_num) {
                case 0:  d->uncommitted  = arg; break;
                case 1:  d->dst_template = arg; break;
                default: argp_usage(state);
            }
            break;

        case ARGP_KEY_END:
            if (state->arg_num < 1) {
                argp_usage(state);
            }
            if (state->arg_num == 1 && d->retain > 1) {
                argp_error(state,
                           "Cannot specify --retain without DESTTEMPLATE");
            }
            if (d->nfields <= 0 && d->nplanes <= 0 && d->nlines <= 0) {
                argp_error(state,
                           "At least one field, plane, or line problem size"
                           " must be specified");
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

        case 'R':
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%d %c", &d->retain, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "retain option is malformed: '%s'", arg);
            }
            if (d->retain < 1) {
                argp_failure(state, EX_USAGE, 0,
                        "retain value %d must be strictly positive",
                        d->retain);
            }
            break;

        case NFIELDS:
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%d %c", &d->nfields, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "nfields is malformed: '%s'", arg);
            }
            if (d->nfields < 0) {
                argp_failure(state, EX_USAGE, 0,
                        "nfields value %d must be nonnegative", d->nfields);
            }
            break;

        case NPLANES:
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%d %c", &d->nplanes, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "nplanes is malformed: '%s'", arg);
            }
            if (d->nplanes < 0) {
                argp_failure(state, EX_USAGE, 0,
                        "nplanes value %d must be nonnegative", d->nplanes);
            }
            break;

        case NLINES:
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%d %c", &d->nlines, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "nlines is malformed: '%s'", arg);
            }
            if (d->nlines < 0) {
                argp_failure(state, EX_USAGE, 0,
                        "nlines value %d must be nonnegative", d->nlines);
            }
            break;

        case FIELD_LAYOUT:
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%d %c", &d->f->layout, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "field-layout option is malformed: '%s'", arg);
            }
            if (d->f->layout < 0 || d->f->layout > esio_field_layout_count()) {
                argp_failure(state, EX_USAGE, 0,
                        "field-layout value %d must be in range [%d, %d]",
                        d->f->layout, 0, esio_field_layout_count());
            }
            break;

        case 'd':
            d->typesize = sizeof(double);
            break;

        case 's':
            d->typesize = sizeof(float);
            break;

        case FIELD_NCOMPONENTS:
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%d %c",
                            &d->f->ncomponents, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "field-ncomponents option is malformed: '%s'", arg);
            }
            if (d->f->ncomponents < 1) {
                argp_failure(state, EX_USAGE, 0,
                        "field-ncomponents value %d must be strictly positive",
                        d->f->ncomponents);
            }
            break;

        case PLANE_NCOMPONENTS:
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%d %c",
                            &d->p->ncomponents, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "plane-ncomponents option is malformed: '%s'", arg);
            }
            if (d->p->ncomponents < 1) {
                argp_failure(state, EX_USAGE, 0,
                        "plane-ncomponents value %d must be strictly positive",
                        d->p->ncomponents);
            }
            break;

        case LINE_NCOMPONENTS:
            errno = 0;
            if (1 != sscanf(arg ? arg : "", "%d %c",
                            &d->l->ncomponents, &ignore)) {
                argp_failure(state, EX_USAGE, errno,
                        "line-ncomponents option is malformed: '%s'", arg);
            }
            if (d->l->ncomponents < 1) {
                argp_failure(state, EX_USAGE, 0,
                        "line-ncomponents value %d must be strictly positive",
                        d->l->ncomponents);
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
            d->nfields = max(d->nfields, 1); // Set nfields >= 1
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
            d->nplanes = max(d->nplanes, 1); // Set nplanes >= 1
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
            d->nlines = max(d->nlines, 1); // Set nlines >= 1
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
            d->nfields = max(d->nfields, 1); // Set nfields >= 1
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
            d->nplanes = max(d->nplanes, 1); // Set nplanes >= 1
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
            d->nlines = max(d->nlines, 1); // Set nlines >= 1
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
    GRVY_TIMER_INIT(argp_program_version);

    // Initialize default argument storage and default values
    struct details d;        memset(&d, 0, sizeof(struct details));
    struct field_details f;  memset(&f, 0, sizeof(struct field_details));
    struct plane_details p;  memset(&p, 0, sizeof(struct plane_details));
    struct line_details  l;  memset(&l, 0, sizeof(struct line_details));
    d.typesize = sizeof(double);
    d.repeat = 1;
    d.retain = 1;
    f.ncomponents = p.ncomponents = l.ncomponents = 1;
    d.f = &f;
    d.p = &p;
    d.l = &l;

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

    // Print program banner that shows version and program arguments
    fprintf(rankout, "%s invoked as\n\t", argp_program_version);
    for (int i = 0; i < argc; ++i) {
        fprintf(rankout, " %s", argv[i]);
    }
    fprintf(rankout, "\n");

    // Initialize ESIO handle
    d.h = esio_handle_initialize(MPI_COMM_WORLD);

    // Initialize the field, plane, and line problems
    if (d.nfields) field_initialize(&d, &f);
    if (d.nplanes) plane_initialize(&d, &p);
    if (d.nlines)  line_initialize( &d, &l);

    // Prepare scratch space to hold names for fields, planes, and lines
    char  *names   = NULL;
    size_t namelen = 0;    // Includes null terminator
    {
        int m = d.nfields;
        m     = (m > d.nplanes) ? m : d.nplanes;
        m     = (m > d.nlines)  ? m : d.nlines;
        int ndigits = ceil(log(m)/log(10));
        namelen = 1 /* f|p|l */ + ndigits + 1 /* null */;
        names = malloc(m * namelen);
        assert(names);
        char *n = names;
        for (int i = 0; i < m; ++i) {
            sprintf(n, "X%0*d", ndigits, i);
            n += namelen;
        }
    }

    // Determine combined size of problem across all ranks
    long localbytes = f.bytes + p.bytes + l.bytes;
    long globalbytes;
    ESIO_MPICHKQ(MPI_Allreduce(&localbytes, &globalbytes, 1,
                               MPI_LONG, MPI_SUM, MPI_COMM_WORLD));
    {
        double coeff;
        const char *units;
        to_human_readable_byte_count(globalbytes, 0, &coeff, &units);
        fprintf(rankout, "Global overall problem size is %.3f %s\n",
                coeff, units);
    }

    fprintf(rankout, "Allocating and filling required memory buffers\n");
    if (d.nfields) f.data = malloc_and_fill(&d, f.bytes);
    if (d.nplanes) p.data = malloc_and_fill(&d, p.bytes);
    if (d.nlines)  l.data = malloc_and_fill(&d, l.bytes);

    // Determine which functions to invoke below based on d.typesize
    int (*p_esio_field_writev)(const esio_handle, const char *,
                               const void *, int, int, int, int) = NULL;
    int (*p_esio_plane_writev)(const esio_handle, const char *,
                               const void *, int, int, int) = NULL;
    int (*p_esio_line_writev)(const esio_handle, const char *,
                              const void *, int, int) = NULL;
    // TODO Understand "warning: assignment from incompatible pointer type"
    // which arise from the function pointer assignments that follow.
    switch (d.typesize)
    {
        case sizeof(double):
            p_esio_field_writev = &esio_field_writev_double;
            p_esio_plane_writev = &esio_plane_writev_double;
            p_esio_line_writev  = &esio_line_writev_double;
            break;
        case sizeof(float):
            p_esio_field_writev = &esio_field_writev_float;
            p_esio_plane_writev = &esio_plane_writev_float;
            p_esio_line_writev  = &esio_line_writev_float;
            break;
        default:
            MPI_Abort(MPI_COMM_WORLD, 1); // Sanity failure
    }

    fprintf(rankout, "Beginning benchmark...\n");
    ESIO_MPICHKQ(MPI_Barrier(MPI_COMM_WORLD)); // Synchronize
    const double start = MPI_Wtime();

    char *n;
    GRVY_TIMER_RESET();
    for (int i = 0; i < d.repeat; ++i) {
        fprintf(rankout, "\tIteration %d\n", i);

        GRVY_TIMER_BEGIN("esio_file_create");
        esio_file_create(d.h, d.uncommitted, 1 /*overwrite*/);
        GRVY_TIMER_END("esio_file_create");

        if (d.nfields) {
            GRVY_TIMER_BEGIN("esio_field_write");
            n = names;
            for (int j = 0; j < d.nfields; ++j) {
                *n = 'f';
                p_esio_field_writev(d.h, n,
                        f.data + j * f.stride, 0, 0, 0, f.ncomponents);
                n += namelen;
            }
            GRVY_TIMER_END("esio_field_write");
        }

        if (d.nplanes) {
            GRVY_TIMER_BEGIN("esio_plane_write");
            n = names;
            for (int j = 0; j < d.nplanes; ++j) {
                *n = 'p';
                p_esio_plane_writev(d.h, n,
                        p.data + j * p.stride, 0, 0, p.ncomponents);
                n += namelen;
            }
            GRVY_TIMER_END("esio_plane_write");
        }

        if (d.nlines) {
            GRVY_TIMER_BEGIN("esio_line_write");
            n = names;
            for (int j = 0; j < d.nlines; ++j) {
                *n = 'l';
                p_esio_line_writev(d.h, n,
                        l.data + j * l.stride, 0, l.ncomponents);
                n += namelen;
            }
            GRVY_TIMER_END("esio_line_write");
        }

        GRVY_TIMER_BEGIN("esio_file_flush");
        esio_file_flush(d.h);
        GRVY_TIMER_END("esio_file_flush");

        if (d.dst_template) {
            GRVY_TIMER_BEGIN("esio_file_close_restart");
            esio_file_close_restart(d.h, d.dst_template, d.retain);
            GRVY_TIMER_END("esio_file_close_restart");
        } else {
            GRVY_TIMER_BEGIN("esio_file_close");
            esio_file_close(d.h);
            GRVY_TIMER_END("esio_file_close");
        }
    }
    GRVY_TIMER_FINALIZE();

    const double end = MPI_Wtime();
    ESIO_MPICHKQ(MPI_Barrier(MPI_COMM_WORLD)); // Synchronize
    fprintf(rankout, "Ending benchmark...\n");

    const double elapsed = end - start;
    const double mean = elapsed / d.repeat;
    {
        double coeff;
        const char *units;
        to_human_readable_byte_count(
                floor(globalbytes / mean), 0, &coeff, &units);
        fprintf(rankout,
            "Mean global transfer rate across %d iteration(s) was %.4f %s/s\n",
            d.repeat, coeff, units);
    }

    // TODO Get timing information back from multiple ranks
    if (d.world_rank == 0) GRVY_TIMER_SUMMARIZE();

    // Finalize the field, plane, and line problems
    if (d.nfields) field_finalize(&d, &f);
    if (d.nplanes) plane_finalize(&d, &p);
    if (d.nlines)  line_finalize( &d, &l);

    // Finalize ESIO handle
    esio_handle_finalize(d.h);

    return 0;
}

void print_version(FILE *stream, struct argp_state *state)
{
    (void) state; // Unused

    fputs(argp_program_version, stream);
    unsigned majnum, minnum, relnum;
    if (H5get_libversion(&majnum, &minnum, &relnum) >= 0) {
        fprintf(stream, " linked against HDF5 %u.%u.%u",
                         majnum, minnum, relnum);
    }
    int version, subversion;
    if (MPI_SUCCESS == MPI_Get_version(&version, &subversion)) {
        fprintf(stream, " running atop MPI %d.%d", version, subversion);
    }
    fputc('\n', stream);
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


static inline int max(int a, int b)
{
    return a > b ? a : b;
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
    int offset = 0;
    for (int i = 0; i < rank; ++i) {
        offset += local(nglobal, nranks, rank, extralow);
    }
    return offset;
}


static int global_minmax(long val, long *min, long *max)
{
    long send[2] = { -val, val };
    long recv[2];
    ESIO_MPICHKQ(MPI_Allreduce(
                send, recv, 2, MPI_LONG, MPI_MAX, MPI_COMM_WORLD));
    *min = -recv[0];
    *max =  recv[1];

    return ESIO_SUCCESS;
}


static int field_initialize(struct details *d, struct field_details *f)
{
    double coeff;
    const char *units;

    fprintf(rankout, "Initializing field problem...\n");

    // Use MPI (temporarily) to find topology for field problem
    ESIO_MPICHKQ(MPI_Dims_create(d->world_size, 3, f->dims));
    MPI_Comm tmp;
    int periods[3] = { 0, 0, 0 };
    ESIO_MPICHKQ(MPI_Cart_create(
                MPI_COMM_WORLD, 3, f->dims, periods, 0, &tmp));
    ESIO_MPICHKQ(MPI_Comm_rank(tmp, &f->rank));
    ESIO_MPICHKQ(MPI_Cart_coords(tmp, f->rank, 3, f->coords));
    ESIO_MPICHKQ(MPI_Comm_free(&tmp));

    fprintf(rankout,
            "\tField spread across (%d x %d x %d)-rank MPI topology\n",
            f->dims[0], f->dims[1], f->dims[2]);

    // Compute global problem size, if necessary, from memory constraint
    if (f->bytes) {
        const double nvectors = (f->bytes * d->world_size)
                              / ((double) f->ncomponents * d->typesize)
                              / ((double) d->nfields);
        f->cglobal = f->bglobal = f->aglobal = ceil(cbrt(nvectors));

        to_human_readable_byte_count(f->bytes, 0, &coeff, &units);
        fprintf(rankout,
            "\tPer-rank %.2f %s memory requested => %d x %d x %d problem\n",
            coeff, units, f->cglobal, f->bglobal, f->aglobal);
    }

    // Ensure a non-trivial problem was requested
    if (!f->cglobal || !f->bglobal || !f->aglobal) {
        fprintf(rankout,
                "You must specify a non-trivial field when nfields > 0!\n");
        MPI_Abort(MPI_COMM_WORLD, EX_USAGE);
    }

    // Compute local portion of global problem
    f->clocal = local(f->cglobal, f->dims[0], f->coords[0], 1);
    f->cstart = start(f->cglobal, f->dims[0], f->coords[0], 1);
    f->blocal = local(f->bglobal, f->dims[1], f->coords[1], 0);
    f->bstart = start(f->bglobal, f->dims[1], f->coords[1], 0);
    f->alocal = local(f->aglobal, f->dims[2], f->coords[2], 1);
    f->astart = start(f->aglobal, f->dims[2], f->coords[2], 1);

    // Compute memory requirement for local portion
    f->stride = f->clocal * f->blocal * f->alocal
              * f->ncomponents * d->typesize;
    f->bytes  = f->stride * d->nfields;

    // Output minimum and maximum memory required across all ranks
    long minbytes, maxbytes;
    global_minmax(f->bytes, &minbytes, &maxbytes); // Allreduce
    to_human_readable_byte_count(minbytes, 0, &coeff, &units);
    fprintf(rankout, "\tMinimum per-rank field memory is %.3f %s\n",
            coeff, units);
    to_human_readable_byte_count(maxbytes, 0, &coeff, &units);
    fprintf(rankout, "\tMaximum per-rank field memory is %.3f %s\n",
            coeff, units);

    // Establish field problem decomposition within ESIO handle
    fprintf(rankout,
            "\tEstablishing field problem within ESIO using layout %d\n",
            f->layout);
    esio_field_layout_set(d->h, f->layout);
    esio_field_establish(d->h, f->cglobal, f->cstart, f->clocal,
                               f->bglobal, f->bstart, f->blocal,
                               f->aglobal, f->astart, f->alocal);

    fprintf(rankout, "Problem contains %d field(s) of size %d x %d x %d"
            " each with %d %zu-byte component(s)\n", d->nfields,
            f->cglobal, f->bglobal, f->aglobal, f->ncomponents,
            d->typesize);

    return ESIO_SUCCESS;
}


static int plane_initialize(struct details *d, struct plane_details *p)
{
    double coeff;
    const char *units;

    fprintf(rankout, "Initializing plane problem...\n");

    // Use MPI (temporarily) to find topology for field problem
    ESIO_MPICHKQ(MPI_Dims_create(d->world_size, 2, p->dims));
    MPI_Comm tmp;
    int periods[2] = { 0, 0 };
    ESIO_MPICHKQ(MPI_Cart_create(
                MPI_COMM_WORLD, 2, p->dims, periods, 0, &tmp));
    ESIO_MPICHKQ(MPI_Comm_rank(tmp, &p->rank));
    ESIO_MPICHKQ(MPI_Cart_coords(tmp, p->rank, 2, p->coords));
    ESIO_MPICHKQ(MPI_Comm_free(&tmp));

    fprintf(rankout, "\tPlane spread across (%d x %d)-rank MPI topology\n",
            p->dims[0], p->dims[1]);

    // Compute global problem size, if necessary, from memory constraint
    if (p->bytes) {
        const double nvectors = (p->bytes * d->world_size)
                              / ((double) p->ncomponents * d->typesize)
                              / ((double) d->nplanes);
        p->bglobal = p->aglobal = ceil(sqrt(nvectors));

        to_human_readable_byte_count(p->bytes, 0, &coeff, &units);
        fprintf(rankout,
            "\tPer-rank %.2f %s memory requested => %d x %d problem\n",
            coeff, units, p->bglobal, p->aglobal);
    }

    // Ensure a non-trivial problem was requested
    if (!p->bglobal || !p->aglobal) {
        fprintf(rankout,
                "You must specify a non-trivial plane when nplanes > 0!\n");
        MPI_Abort(MPI_COMM_WORLD, EX_USAGE);
    }

    // Compute local portion of global problem
    p->blocal = local(p->bglobal, p->dims[0], p->coords[0], 0);
    p->bstart = start(p->bglobal, p->dims[0], p->coords[0], 0);
    p->alocal = local(p->aglobal, p->dims[1], p->coords[1], 1);
    p->astart = start(p->aglobal, p->dims[1], p->coords[1], 1);

    // Compute memory requirement for local portion
    p->stride = p->blocal * p->alocal * p->ncomponents * d->typesize;
    p->bytes  = p->stride * d->nplanes;

    // Output minimum and maximum memory required across all ranks
    long minbytes, maxbytes;
    global_minmax(p->bytes, &minbytes, &maxbytes); // Allreduce
    to_human_readable_byte_count(minbytes, 0, &coeff, &units);
    fprintf(rankout, "\tMinimum per-rank plane memory is %.3f %s\n",
            coeff, units);
    to_human_readable_byte_count(maxbytes, 0, &coeff, &units);
    fprintf(rankout, "\tMaximum per-rank plane memory is %.3f %s\n",
            coeff, units);

    // Establish plane problem decomposition within ESIO handle
    fprintf(rankout, "\tEstablishing plane problem within ESIO\n");
    esio_plane_establish(d->h, p->bglobal, p->bstart, p->blocal,
                               p->aglobal, p->astart, p->alocal);

    fprintf(rankout, "Problem contains %d plane(s) of size %d x %d"
            " each with %d %zu-byte component(s)\n", d->nplanes,
            p->bglobal, p->aglobal, p->ncomponents, d->typesize);

    return ESIO_SUCCESS;
}


static int line_initialize(struct details *d, struct line_details  *l)
{
    double coeff;
    const char *units;

    fprintf(rankout, "Initializing line problem...\n");

    // Use MPI (temporarily) to find topology for field problem
    ESIO_MPICHKQ(MPI_Dims_create(d->world_size, 1, l->dims));
    MPI_Comm tmp;
    int periods[1] = { 0 };
    ESIO_MPICHKQ(MPI_Cart_create(
                MPI_COMM_WORLD, 1, l->dims, periods, 0, &tmp));
    ESIO_MPICHKQ(MPI_Comm_rank(tmp, &l->rank));
    ESIO_MPICHKQ(MPI_Cart_coords(tmp, l->rank, 1, l->coords));
    ESIO_MPICHKQ(MPI_Comm_free(&tmp));

    fprintf(rankout, "\tLine spread across %d-rank MPI topology\n",
            l->dims[0]);

    // Compute global problem size, if necessary, from memory constraint
    if (l->bytes) {
        const double nvectors = (l->bytes * d->world_size)
                              / ((double) l->ncomponents * d->typesize)
                              / ((double) d->nlines);
        l->aglobal = ceil(nvectors);

        to_human_readable_byte_count(l->bytes, 0, &coeff, &units);
        fprintf(rankout,
            "\tPer-rank %.2f %s memory requested => %d problem\n",
            coeff, units, l->aglobal);
    }

    // Ensure a non-trivial problem was requested
    if (!l->aglobal) {
        fprintf(rankout,
                "You must specify a non-trivial line when nlines > 0!\n");
        MPI_Abort(MPI_COMM_WORLD, EX_USAGE);
    }

    // Compute local portion of global problem
    l->alocal = local(l->aglobal, l->dims[0], l->coords[0], 1);
    l->astart = start(l->aglobal, l->dims[0], l->coords[0], 1);

    // Compute memory requirement for local portion
    l->stride = l->alocal * l->ncomponents * d->typesize;
    l->bytes  = l->stride * d->nlines;

    // Output minimum and maximum memory required across all ranks
    long minbytes, maxbytes;
    global_minmax(l->bytes, &minbytes, &maxbytes); // Allreduce
    to_human_readable_byte_count(minbytes, 0, &coeff, &units);
    fprintf(rankout, "\tMinimum per-rank line memory is %.3f %s\n",
            coeff, units);
    to_human_readable_byte_count(maxbytes, 0, &coeff, &units);
    fprintf(rankout, "\tMaximum per-rank line memory is %.3f %s\n",
            coeff, units);

    // Establish line problem decomposition within ESIO handle
    fprintf(rankout, "\tEstablishing line problem within ESIO\n");
    esio_line_establish(d->h, l->aglobal, l->astart, l->alocal);

    fprintf(rankout, "Problem contains %d line(s) of size %d"
            " each with %d %zu-byte component(s)\n", d->nlines,
            l->aglobal, l->ncomponents, d->typesize);

    return ESIO_SUCCESS;
}


static void* malloc_and_fill(struct details *d, const long bytes)
{
    // Malloc
    void *p = malloc(bytes);
    if (!p) {
        double coeff;
        const char *units;
        to_human_readable_byte_count(bytes, 0, &coeff, &units);
        fprintf(stderr, "Unable to malloc %.2f %s bytes on rank %d\n",
                coeff, units, d->world_rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Store previous random seed and provide a new one
    char state[64];
    initstate(d->world_rank, state, sizeof(state)/sizeof(state[0]));
    char *previous = setstate(state);

    // Fill
    const size_t count = bytes / d->typesize;
    switch (d->typesize)
    {
        case sizeof(double):
            for (size_t i = 0; i < count; ++i)
                ((double *) p)[i] = (double) random();
            break;
        case sizeof(float):
            for (size_t i = 0; i < count; ++i)
                ((float *) p)[i] = (float) random();
            break;
    }

    // Restore previous random state
    setstate(previous);

    return p;
}


static int field_finalize(struct details *d, struct field_details *f)
{
    (void) d; // Unused

    if (f && f->data) {
        free(f->data);
        f->data = NULL;
    }

    return ESIO_SUCCESS;
}


static int plane_finalize(struct details *d, struct plane_details *p)
{
    (void) d; // Unused

    if (p && p->data) {
        free(p->data);
        p->data = NULL;
    }

    return ESIO_SUCCESS;
}


static int line_finalize( struct details *d, struct line_details  *l)
{
    (void) d; // Unused

    if (l && l->data) {
        free(l->data);
        l->data = NULL;
    }

    return ESIO_SUCCESS;
}
