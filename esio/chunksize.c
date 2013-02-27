//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// ESIO 0.1.9: ExaScale IO library for turbulence simulation restart files
// http://red.ices.utexas.edu/projects/esio/
//
// Copyright (C) 2010, 2011, 2012, 2013 The PECOS Development Team
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
#include "chunksize.h"

#include <mpi.h>

#include "error.h"

int
chunksize_line(MPI_Comm comm,
               int aglobal, int astart, int alocal, int *achunk)
{
    // Ignored parameters kept for API consistency
    (void) astart;

    // Obtain maximum and minimum aglobal, alocal values in a single Allreduce
    int sendbuf[] = { aglobal, -aglobal, alocal, -alocal };
    const int count = sizeof(sendbuf)/sizeof(sendbuf[0]);
    int recvbuf[count];
    ESIO_MPICHKQ(MPI_Allreduce(sendbuf, recvbuf, count,
                               MPI_INT, MPI_MAX, comm));
    const int aglobal_max =  recvbuf[0];
    const int aglobal_min = -recvbuf[1];
    const int alocal_max  =  recvbuf[2];
    const int alocal_min  = -recvbuf[3];

    // Check that global values match on all ranks.
    // Hides usage error checking costs in other global communication work.
    if (aglobal_max != aglobal_min) {
        ESIO_ERROR("Value supplied for aglobal varies by rank", ESIO_EINVAL);
    }

    // Compute the chunk size to use based on global information
    (void) alocal_min; // For now the minimum value is ignored
    *achunk = alocal_max;

    return ESIO_SUCCESS;
}

int
chunksize_plane(MPI_Comm comm,
                int bglobal, int bstart, int blocal, int *bchunk,
                int aglobal, int astart, int alocal, int *achunk)
{
    // Ignored parameters kept for API consistency
    (void) bstart;
    (void) astart;

    // Obtain maximum and minimum global, local values in a single Allreduce
    int sendbuf[] = {
        bglobal, -bglobal, blocal, -blocal,
        aglobal, -aglobal, alocal, -alocal
    };
    const int count = sizeof(sendbuf)/sizeof(sendbuf[0]);
    int recvbuf[count];
    ESIO_MPICHKQ(MPI_Allreduce(sendbuf, recvbuf, count,
                               MPI_INT, MPI_MAX, comm));
    const int bglobal_max =  recvbuf[0];
    const int bglobal_min = -recvbuf[1];
    const int blocal_max  =  recvbuf[2];
    const int blocal_min  = -recvbuf[3];
    const int aglobal_max =  recvbuf[4];
    const int aglobal_min = -recvbuf[5];
    const int alocal_max  =  recvbuf[6];
    const int alocal_min  = -recvbuf[7];

    // Check that global values match on all ranks.
    // Hides usage error checking costs in other global communication work.
    if (bglobal_max != bglobal_min) {
        ESIO_ERROR("Value supplied for bglobal varies by rank", ESIO_EINVAL);
    }
    if (aglobal_max != aglobal_min) {
        ESIO_ERROR("Value supplied for aglobal varies by rank", ESIO_EINVAL);
    }

    // Compute the chunk size to use based on global information
    (void) blocal_min; // For now the minimum value is ignored
    *bchunk = blocal_max;
    (void) alocal_min; // For now the minimum value is ignored
    *achunk = alocal_max;

    return ESIO_SUCCESS;
}

int
chunksize_field(MPI_Comm comm,
                int cglobal, int cstart, int clocal, int *cchunk,
                int bglobal, int bstart, int blocal, int *bchunk,
                int aglobal, int astart, int alocal, int *achunk)
{
    // Ignored parameters kept for API consistency
    (void) cstart;
    (void) bstart;
    (void) astart;

    // Obtain maximum and minimum global, local values in a single Allreduce
    int sendbuf[] = {
        cglobal, -cglobal, clocal, -clocal,
        bglobal, -bglobal, blocal, -blocal,
        aglobal, -aglobal, alocal, -alocal
    };
    const int count = sizeof(sendbuf)/sizeof(sendbuf[0]);
    int recvbuf[count];
    ESIO_MPICHKQ(MPI_Allreduce(sendbuf, recvbuf, count,
                               MPI_INT, MPI_MAX, comm));
    const int cglobal_max =  recvbuf[ 0];
    const int cglobal_min = -recvbuf[ 1];
    const int clocal_max  =  recvbuf[ 2];
    const int clocal_min  = -recvbuf[ 3];
    const int bglobal_max =  recvbuf[ 4];
    const int bglobal_min = -recvbuf[ 5];
    const int blocal_max  =  recvbuf[ 6];
    const int blocal_min  = -recvbuf[ 7];
    const int aglobal_max =  recvbuf[ 8];
    const int aglobal_min = -recvbuf[ 9];
    const int alocal_max  =  recvbuf[10];
    const int alocal_min  = -recvbuf[11];

    // Check that global values match on all ranks.
    // Hides usage error checking costs in other global communication work.
    if (cglobal_max != cglobal_min) {
        ESIO_ERROR("Value supplied for cglobal varies by rank", ESIO_EINVAL);
    }
    if (bglobal_max != bglobal_min) {
        ESIO_ERROR("Value supplied for bglobal varies by rank", ESIO_EINVAL);
    }
    if (aglobal_max != aglobal_min) {
        ESIO_ERROR("Value supplied for aglobal varies by rank", ESIO_EINVAL);
    }

    // Compute the chunk size to use based on global information
    (void) clocal_min; // For now the minimum value is ignored
    *cchunk = clocal_max;
    (void) blocal_min; // For now the minimum value is ignored
    *bchunk = blocal_max;
    (void) alocal_min; // For now the minimum value is ignored
    *achunk = alocal_max;

    return ESIO_SUCCESS;
}
