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

program line_int_f

    use mpi
    use esio
    use testframework

    implicit none

    integer :: i, j
    integer :: value_scalar(3),   buffer_scalar(3)
    integer :: value_vector(2,3), buffer_vector(2,3)
    integer :: aglobal, astart, alocal, astride, ncomponents

    call testframework_setup()

!   Populate test data based on process rank
    aglobal     = size(value_scalar) * world_size
    astart      = size(value_scalar) * world_rank + 1
    alocal      = size(value_scalar)
    astride     = 1
    ncomponents = size(value_vector,1)
    do j = lbound(value_scalar,1), ubound(value_scalar,1)
       value_scalar(j) = j + astart - 1
    end do
    do j = lbound(value_vector,2), ubound(value_vector,2)
        do i = lbound(value_vector,1), ubound(value_vector,1)
            value_vector(i,j) = i*value_scalar(j)
        end do
    end do

!   Create a file
    call esio_file_create(h, filename, .false., ierr)
    ASSERT(ierr == 0)

!   Write a scalar-valued line
    call esio_line_write_integer(h, "name_scalar", value_scalar,   &
                                 aglobal, astart, alocal, astride, &
                                 ierr)
    ASSERT(ierr == 0)

!   Write a vector-valued line
    call esio_line_writev_integer(h, "name_vector", value_vector, &
                                  ncomponents,                    &
                                  aglobal, astart, alocal, 0,     &
                                  ierr)
    ASSERT(ierr == 0)

!   Close the file
    call esio_file_close(h, ierr)
    ASSERT(ierr == 0)

!   Re-open the file read-only
    call esio_file_open(h, filename, .false., ierr)
    ASSERT(ierr == 0)

!   Check the scalar-valued line's size and data
    call esio_line_size(h, "name_scalar", i, ierr)
    ASSERT(ierr == 0)
    ASSERT(i == size(value_scalar))
    call esio_line_sizev(h, "name_scalar", i, j, ierr)
    ASSERT(ierr == 0)
    ASSERT(i == 1)
    ASSERT(j == size(value_scalar))
    call esio_line_read_integer(h, "name_scalar", buffer_scalar,  &
                                aglobal, astart, alocal, astride, &
                                ierr)
    ASSERT(ierr == 0)
    ASSERT(all(buffer_scalar == value_scalar))

!   Check the vector-valued line's size and data
    call esio_line_sizev(h, "name_vector", i, j, ierr)
    ASSERT(ierr == 0)
    ASSERT(i == ncomponents)
    ASSERT(j == aglobal)
    call esio_line_readv_integer(h, "name_vector", buffer_vector,  &
                                 ncomponents,                      &
                                 aglobal, astart, alocal, 0,       &
                                 ierr)
    ASSERT(ierr == 0)
    ASSERT(all(buffer_vector == value_vector))

!   Close the file
    call esio_file_close(h, ierr)
    ASSERT(ierr == 0)

    call testframework_teardown()

end program line_int_f
