//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// esio 0.1.5: ExaScale IO library for turbulence simulation restart files
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

#ifndef __ESIO_MPI_ARGP_H__
#define __ESIO_MPI_ARGP_H__

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

#endif /* __ESIO_MPI_ARGP_H__ */
