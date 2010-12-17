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
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <esio/restart-rename.h>

static void print_usage(FILE *stream, const char *arg0)
{
    fprintf(stream, "Usage:\n");
    fprintf(stream, "\t%s src_filename dst_template retain_count\n", arg0);
}

int main(int argc, char *argv[])
{
    int c;
    while ((c = getopt(argc, argv, "h")) != -1) {
        switch (c) {
            case 'h':
                print_usage(stdout, argv[0]);
                return EXIT_SUCCESS;
            case '?':
                if (isprint(optopt)) {
                    fprintf(stderr, "Unknown option '-%c'.\n", c);
                }
                print_usage(stderr, argv[0]);
                return EXIT_FAILURE;

        }
    }
    if (argc < 4) {
        print_usage(stderr, argv[0]);
        return EXIT_FAILURE;
    }

    const char *src_filename = argv[1];
    const char *dst_template = argv[2];
    errno = 0;
    const long nlong = strtol(argv[3], NULL, 10);
    if (errno) {
        fprintf(stderr, "Error parsing retain_count = '%s'\n", argv[3]);
        return EXIT_FAILURE;
    }
    if (nlong < 1) {
        fprintf(stderr, "Error retain_count < 1\n");
        return EXIT_FAILURE;
    }
    if (nlong > INT_MAX) {
        fprintf(stderr, "Error parsing retain_count > INT_MAX\n");
        return EXIT_FAILURE;
    }
    const int n = (int) nlong;

    return restart_rename(src_filename, dst_template, n);
}
