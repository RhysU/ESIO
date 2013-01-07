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

#ifndef ESIO_MPI_ARGP_H__
#define ESIO_MPI_ARGP_H__

#include "argp.h"

/**
 * Call <a href="http://www.gnu.org/s/libc/manual/html_node/Argp.html">
 * Argp</a>'s \c argp_parse in an MPI-friendly way.  Processes with
 * nonzero rank will have their \c stdout and \c stderr redirected
 * to <tt>/dev/null</tt> during \c argp_parse.
 *
 * @param rank MPI rank of this process.  Output from \c argp_parse
 *             will only be observable from rank zero.
 * @param argp      Per \c argp_parse semantics.
 * @param argc      Per \c argp_parse semantics.
 * @param argv      Per \c argp_parse semantics.
 * @param flags     Per \c argp_parse semantics.
 * @param arg_index Per \c argp_parse semantics.
 * @param input     Per \c argp_parse semantics.
 *
 * @return Per \c argp_parse semantics.
 */
error_t mpi_argp_parse(const int rank,
                       const struct argp *argp,
                       int argc,
                       char **argv,
                       unsigned flags,
                       int *arg_index,
                       void *input);

#endif /* ESIO_MPI_ARGP_H__ */
