//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// esio 0.1.3: ExaScale IO library for turbulence simulation restart files
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

#include <stdio.h>
#include <unistd.h>
#include "argp.h"
#include "mpi_argp.h"

error_t mpi_argp_parse(const int rank,
                       const struct argp *argp,
                       int argc,
                       char **argv,
                       unsigned flags,
                       int *arg_index,
                       void *input)
{
    // Flush stdout, stderr
    if (fflush(stdout))
        perror("mpi_argp_parse error flushing stdout prior to redirect");
    if (fflush(stderr))
        perror("mpi_argp_parse Error flushing stderr prior to redirect");

    // Save stdout, stderr so we may restore them later
    int stdout_copy, stderr_copy;
    if ((stdout_copy = dup(fileno(stdout))) < 0)
        perror("mpi_argp_parse error duplicating stdout");
    if ((stderr_copy = dup(fileno(stderr))) < 0)
        perror("mpi_argp_parse error duplicating stderr");

    // On non-root processes redirect stdout, stderr to /dev/null
    if (rank) {
        if (!freopen("/dev/null", "a", stdout))
            perror("mpi_argp_parse error redirecting stdout");
        if (!freopen("/dev/null", "a", stderr))
            perror("mpi_argp_parse error redirecting stderr");
    }

    // Invoke argp per http://www.gnu.org/s/libc/manual/html_node/Argp.html
    error_t retval = argp_parse(argp, argc, argv, flags, arg_index, input);

    // Flush stdout, stderr again
    if (fflush(stdout))
        perror("mpi_argp_parse error flushing stdout after redirect");
    if (fflush(stderr))
        perror("mpi_argp_parse error flushing stderr after redirect");

    // Restore stdout, stderr
    if (dup2(stdout_copy, fileno(stdout)) < 0)
        perror("mpi_argp_parse error reopening stdout");
    if (dup2(stderr_copy, fileno(stderr)) < 0)
        perror("mpi_argp_parse error reopening stderr");

    // Close saved versions of stdout, stderr
    if (close(stdout_copy))
        perror("mpi_argp_parse error closing stdout_copy");
    if (close(stderr_copy))
        perror("mpi_argp_parse error closing stderr_copy");

    // Clear any errors that may have occurred on stdout, stderr
    clearerr(stdout);
    clearerr(stderr);

    // Return what argp_parse returned
    return retval;
}
