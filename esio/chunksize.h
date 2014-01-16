//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// ESIO 0.2.0: ExaScale IO library for turbulence simulation restart files
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

#ifndef ESIO_CHUNKSIZE_H
#define ESIO_CHUNKSIZE_H

//****************************************************************
// INTERNAL INTERNAL INTERNAL INTERNAL INTERNAL INTERNAL INTERNAL
//****************************************************************

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Collectively deduce and return appropriate HDF5 chunk sizes for
 * a line decomposed according to the given parameters.
 *
 * \param[in]  comm    MPI communicator used for any reduction process.
 * \param[in]  aglobal Global number of scalars within the line.
 * \param[in]  astart  Global starting offset (zero-indexed) handled
 *                     locally by this MPI rank.
 * \param[in]  alocal  Number of scalars this MPI rank will write.
 * \param[out] achunk  Chunk size that should be used.
 *
 * \return Either ESIO_SUCCESS \c (0) or one of ::esio_status on failure.
 *
 * @see <a
 * href="http://www.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-SetChunk">
 * <code>H5Pset_chunk</code></a> for more details on chunking.
 */
int
chunksize_line(MPI_Comm comm,
               int aglobal, int astart, int alocal, int *achunk);

/**
 * Collectively deduce and return appropriate HDF5 chunk sizes for
 * a plane decomposed according to the given parameters.
 *
 * \param[in]  comm    MPI communicator used for any reduction process.
 * \param[in]  bglobal Global number of scalars in the slower "B" direction.
 * \param[in]  bstart  Global starting "B" offset.
 * \param[in]  blocal  Number of scalars in "B" this MPI rank will write.
 * \param[out] bchunk  Chunk size to use in the slower "B" direction.
 * \param[in]  aglobal Global number of scalars in the faster "A" direction.
 * \param[in]  astart  Global starting "A" offset.
 * \param[in]  alocal  Number of scalars in "A" this MPI rank will write.
 * \param[out] achunk  Chunk size to use in the faster "A" direction.
 *
 * \return Either ESIO_SUCCESS \c (0) or one of ::esio_status on failure.
 *
 * @see <a
 * href="http://www.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-SetChunk">
 * <code>H5Pset_chunk</code></a> for more details on chunking.
 */
int
chunksize_plane(MPI_Comm comm,
                int bglobal, int bstart, int blocal, int *bchunk,
                int aglobal, int astart, int alocal, int *achunk);

/**
 * Collectively deduce and return appropriate HDF5 chunk sizes for
 * a field decomposed according to the given parameters.
 *
 * \param[in]  comm MPI communicator used for any reduction process.
 * \param[in]  cglobal Global number of scalars in the "C" slowest direction.
 * \param[in]  cstart  Global starting "C" offset.
 * \param[in]  clocal  Number of scalars in "C" this MPI rank will write.
 * \param[out] cchunk  Chunk size to use in the "A" direction.
 * \param[in]  bglobal Global number of scalars in the "B" direction.
 * \param[in]  bstart  Global starting "B" offset.
 * \param[in]  blocal  Number of scalars in "B" this MPI rank will write.
 * \param[out] bchunk  Chunk size to use in the "A" direction.
 * \param[in]  aglobal Global number of scalars in the fastest "A" direction.
 * \param[in]  astart  Global starting "A" offset.
 * \param[in]  alocal  Number of scalars in "A" this MPI rank will write.
 * \param[out] achunk  Chunk size to use in the "A" direction.
 *
 * \return Either ESIO_SUCCESS \c (0) or one of ::esio_status on failure.
 *
 * @see <a
 * href="http://www.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-SetChunk">
 * <code>H5Pset_chunk</code></a> for more details on chunking.
 */
int
chunksize_field(MPI_Comm comm,
                int cglobal, int cstart, int clocal, int *cchunk,
                int bglobal, int bstart, int blocal, int *bchunk,
                int aglobal, int astart, int alocal, int *achunk);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* ESIO_CHUNKSIZE_H */
