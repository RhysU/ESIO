//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// ESIO 0.2.0: ExaScale IO library for turbulence simulation restart files
// http://github.com/RhysU/ESIO
//
// Copyright (C) 2010-2017 The PECOS Development Team
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

#ifndef ESIO_H5UTILS_H
#define ESIO_H5UTILS_H

#include <hdf5.h>

//****************************************************************
// INTERNAL INTERNAL INTERNAL INTERNAL INTERNAL INTERNAL INTERNAL
//****************************************************************

#ifdef __cplusplus
extern "C" {
#endif

// Macro to save and disable current HDF5 error handler
// "id" is present to allow multiple disable/enable calls in one scope.
#define DISABLE_HDF5_ERROR_HANDLER(id)                                   \
    H5E_auto2_t hdf5_handler##id;                                        \
    void *hdf5_client_data##id;                                          \
    H5Eget_auto2(H5E_DEFAULT, &hdf5_handler##id, &hdf5_client_data##id); \
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

// Macro to restore previously saved HDF5 error handler
// "id" is present to allow multiple disable/enable calls in one scope.
#define ENABLE_HDF5_ERROR_HANDLER(id)                                 \
    H5Eset_auto2(H5E_DEFAULT, hdf5_handler##id, hdf5_client_data##id);

herr_t
esio_H5LTget_attribute_ndims_info(hid_t loc_id,
                                  const char *obj_name,
                                  const char *attr_name,
                                  int *rank,
                                  hsize_t *dims,
                                  H5T_class_t *type_class,
                                  size_t *type_size);

// Walk the current error stack to see if errnum appears anywhere within it
htri_t
esio_H5Equery_stack(hid_t errnum);

#ifdef __cplusplus
}
#endif

#endif /* ESIO_H5UTILS_H */
