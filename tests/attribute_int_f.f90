!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! esio 0.0.1: ExaScale IO library for turbulence simulation restart files
!! http://pecos.ices.utexas.edu/
!!
!! Copyright (C) 2010 The PECOS Development Team
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

    use mpi
    use esio
    use testframework

    implicit none

    integer, parameter :: value_scalar    = 123
    integer, parameter :: value_vector(4) = (/12, 34, 56, 78/)
    integer            :: buffer(4), i

    call testframework_setup()

!   Create a file, write values, and close it
    call esio_file_create(h, filename, .false., ierr)
    ASSERT(ierr == 0)
    call esio_attribute_write_integer(h, "name_scalar", value_scalar, ierr)
    ASSERT(ierr == 0)
    call esio_attribute_writev_integer(h, "name_vector", value_vector, &
                                       size(value_vector), ierr)
    ASSERT(ierr == 0)
    call esio_file_close(h, ierr)
    ASSERT(ierr == 0)

!   Re-open the file read-only, check sizes and values, and close it
    call esio_file_open(h, filename, .false., ierr)
    ASSERT(ierr == 0)
    call esio_attribute_sizev(h, "name_scalar", buffer(1), ierr)
    ASSERT(buffer(1) == 1)
    call esio_attribute_read_integer(h, "name_scalar", buffer(1), ierr)
    ASSERT(ierr == 0)
    ASSERT(buffer(1) == 123)

    call esio_attribute_sizev(h, "name_vector", buffer(1), ierr)
    ASSERT(buffer(1) == size(value_vector))
    call esio_attribute_readv_integer(h, "name_vector", buffer, &
                                      size(value_vector), ierr)
    ASSERT(ierr == 0)
    do i = 1, size(value_vector)
        ASSERT(buffer(i) == value_vector(i))
    end do
    call esio_file_close(h, ierr)
    ASSERT(ierr == 0)

    call testframework_teardown()

end program attribute_int_f
