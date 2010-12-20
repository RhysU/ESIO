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

#include <stdlib.h>

#include "argp.h"

#include <mpi.h>
#include <esio/esio.h>

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    atexit((void (*) ()) MPI_Finalize);
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    esio_handle h = esio_handle_initialize(MPI_COMM_WORLD);
    esio_file_create(h, "data.h5", 1 /* overwrite */);

    esio_string_set(h, "program", argv[0]);

    int version = 1;
    esio_attribute_write_int(h, "version", &version);

    double example[2] = { 2.0 * world_rank, 2.0 * world_rank + 1 };
    esio_line_write_double(h, "example", example,
                           2*world_size, 2*world_rank, 2, 1);

    esio_file_close(h);
    esio_handle_finalize(h);

    return 0;
}
