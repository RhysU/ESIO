//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// ESIO 0.1.9: ExaScale IO library for turbulence simulation restart files
// http://red.ices.utexas.edu/projects/esio/
//
// Copyright (C) 2010-2014 The PECOS Development Team
//
// This file is part of ESIO.
//
// ESIO is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3.0 of the License, or
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>
#include <unistd.h>

#include "argp.h"
#include <esio/restart-rename.h>

// Argp details: http://www.gnu.org/s/libc/manual/html_node/Argp.html
const char *argp_program_version      = "esio_rename " PACKAGE_VERSION;
const char *argp_program_bug_address  = PACKAGE_BUGREPORT;
static const char args_doc[]          = "SOURCE DESTTEMPLATE";
static const char doc[] =
"Rename file SOURCE to match DESTTEMPLATE preserving related restart files"
"\v"
"Rename the restart file SOURCE it to match the path given in "
"DESTTEMPLATE. Up to retain_count previous restart files will be "
"retained and will automatically have their index numbers incremented. "
"Index numbers are in the range [0, retain_count-1] (inclusive) with "
"index 0 being the newest file.\n"
"\n"
"DESTTEMPLATE can be an absolute or relative path where the final filename "
"must contain a sequence of one or more consecutive hash signs ('#') which "
"will be populated with restart index numbers.  A sufficient number of "
"leading zeros to accommodate (retain count - 1) separate restart "
"files will always be present.  Using additional hash signs will increase "
"the number of leading zeros appearing in restart file names.  Note that "
"hash signs need to be escaped when appearing in most shell commands.\n"
"\n"
"Note that SOURCE must not match DESTTEMPLATE otherwise this command will "
"fail with mysterious renaming errors.\n"
;

static struct argp_option options[] = {
    { "retain", 'r', "count", 0, "Number of restart files to retain", 0 },
    { 0,        0,   0,       0,  0,                                  0 }
};

// Used by main to communicate with parse_opt
struct arguments {
    char *src_filename;
    char *dst_template;
    int   retain_count;
};

// Parse a single option following Argp semantics
static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
    // Get the input argument from argp_parse.
    struct arguments *arguments = state->input;

    switch (key) {
        case ARGP_KEY_ARG:
            switch (state->arg_num) {
                case 0:  arguments->src_filename = arg; break;
                case 1:  arguments->dst_template = arg; break;
                default: argp_usage(state);
            }
            break;

        case ARGP_KEY_END:
            if (state->arg_num < 2) argp_usage(state);
            break;

        case 'r':
            {
                errno = 0;
                const long n
                  = arg ? strtol(arg, NULL, 10) : arguments->retain_count;
                if (errno || n < 0 || n > INT_MAX) {
                    argp_failure(state, EX_USAGE, errno,
                                 "retain count must be in range [%d, %d]",
                                 1, INT_MAX - 1);
                } else {
                    arguments->retain_count = (int) n;
                }
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
    struct arguments arguments;
    arguments.src_filename = NULL;
    arguments.dst_template = NULL;
    arguments.retain_count = 10;

    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    return restart_rename(arguments.src_filename,
                          arguments.dst_template,
                          arguments.retain_count);
}
