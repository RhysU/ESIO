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

#ifndef __ESIO_H5UTILS_H
#define __ESIO_H5UTILS_H

//****************************************************************
// INTERNAL INTERNAL INTERNAL INTERNAL INTERNAL INTERNAL INTERNAL
//****************************************************************

// Macro to save and disable current HDF5 error handler
#define DISABLE_HDF5_ERROR_HANDLER                               \
    H5E_auto2_t hdf5_handler;                                    \
    void *hdf5_client_data;                                      \
    H5Eget_auto2(H5E_DEFAULT, &hdf5_handler, &hdf5_client_data); \
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

// Macro to restore previously saved HDF5 error handler
#define ENABLE_HDF5_ERROR_HANDLER                               \
    H5Eset_auto2(H5E_DEFAULT, hdf5_handler, hdf5_client_data);

#endif /* __ESIO_H5UTILS_H */
