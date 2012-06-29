!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! esio 0.1.7: ExaScale IO library for turbulence simulation restart files
!! http://pecos.ices.utexas.edu/
!!
!! Copyright (C) 2010, 2011 The PECOS Development Team
!!
!! This library is free software; you can redistribute it and/or
!! modify it under the terms of the Version 2.1 GNU Lesser General
!! Public License as published by the Free Software Foundation.
!!
!! This library is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!! Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public
!! License along with this library; if not, write to the Free Software
!! Foundation, Inc. 51 Franklin Street, Fifth Floor,
!! Boston, MA  02110-1301  USA
!!
!!-----------------------------------------------------------------------el-
!! $Id$

#include "testframework_assert.h"

program attribute_int_f

    use esio
    use testframework

    implicit none

    integer, parameter :: value_scalar    = 123
    integer, parameter :: value_vector(4) = (/12, 34, 56, 78/)
    integer            :: buffer(4)

    call testframework_setup(__FILE__)

!   Create a file
    call esio_file_create(h, filename, .false., ierr)
    ASSERT(ierr == 0)

!   Write a scalar-valued attribute using an explicitly-typed method
    call esio_attribute_write_integer(  &
            h, "/", "name_scalar", value_scalar, ierr)
    ASSERT(ierr == 0)

!   Overwrite the scalar-valued attribute using a generic interface
    call esio_attribute_write(  &
            h, "/", "name_scalar", value_scalar, ierr)
    ASSERT(ierr == 0)

!   Write a vector-valued attribute using an explicitly-typed method
    call esio_attribute_writev_integer(  &
            h, "/", "name_vector", value_vector, size(value_vector), ierr)
    ASSERT(ierr == 0)

!   Overwrite the vector-valued attribute using a generic interface
    call esio_attribute_writev(  &
            h, "/", "name_vector", value_vector, size(value_vector), ierr)
    ASSERT(ierr == 0)

!   Close the file
    call esio_file_close(h, ierr)
    ASSERT(ierr == 0)

!   Re-open the same file read-only
    call esio_file_open(h, filename, .false., ierr)
    ASSERT(ierr == 0)

!   Check the scalar-valued attribute's size and data
    call esio_attribute_sizev(h, "/", "name_scalar", buffer(1), ierr)
    ASSERT(ierr == 0)
    ASSERT(buffer(1) == 1)
    buffer(:) = 0
    call esio_attribute_read_integer(h, "/", "name_scalar", buffer(1), ierr)
    ASSERT(ierr == 0)
    ASSERT(buffer(1) == 123)
    buffer(:) = 0
    call esio_attribute_read(h, "/", "name_scalar", buffer(1), ierr)
    ASSERT(ierr == 0)
    ASSERT(buffer(1) == 123)
    buffer(:) = 0

!   Check the vector-valued attribute's size and data
    call esio_attribute_sizev(h, "/", "name_vector", buffer(1), ierr)
    ASSERT(ierr == 0)
    ASSERT(buffer(1) == size(value_vector))
    call esio_attribute_sizev(h, "/", "DOESNOTEXIST", buffer(1), ierr)
    ASSERT(ierr /= 0)
    ASSERT(buffer(1) == size(value_vector))  ! Preserved per feature #2419?
    buffer(:) = 0
    call esio_attribute_readv_integer(  &
            h, "/", "name_vector", buffer, size(value_vector), ierr)
    ASSERT(ierr == 0)
    ASSERT(all(buffer == value_vector))
    buffer(:) = 0
    call esio_attribute_readv(  &
            h, "/", "name_vector", buffer, size(value_vector), ierr)
    ASSERT(ierr == 0)
    ASSERT(all(buffer == value_vector))
    buffer(:) = 0

!   Close the file
    call esio_file_close(h, ierr)
    ASSERT(ierr == 0)

    call testframework_teardown()

end program attribute_int_f
