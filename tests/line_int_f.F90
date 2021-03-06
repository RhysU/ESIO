!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! ESIO 0.2.0: ExaScale IO library for turbulence simulation restart files
!! http://github.com/RhysU/ESIO
!!
!! Copyright (C) 2010-2017 The PECOS Development Team
!!
!! This file is part of ESIO.
!!
!! ESIO is free software: you can redistribute it and/or modify
!! it under the terms of the GNU Lesser General Public License as published
!! by the Free Software Foundation, either version 3.0 of the License, or
!! (at your option) any later version.
!!
!! ESIO is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public License
!! along with ESIO.  If not, see <http://www.gnu.org/licenses/>.
!!
!!-----------------------------------------------------------------------el-
!! $Id$

#include "testframework_assert.h"

program line_int_f

    use esio
    use testframework

    implicit none

    integer            :: value_scalar(3),   buffer_scalar(3)
    integer            :: value_vector(2,3), buffer_vector(2,3)
    integer, parameter :: ncomponents = size(value_vector, 1)
    integer            :: i, j

    ! One-dimensional Cartesian topology
    call testframework_setup (__FILE__, 1)

!   Prepare homogeneous test data distribution across the topology
    do i = 1, ndims
        global(i) = size(value_scalar,i) * dims(i)
        start(i)  = size(value_scalar,i) * coords(i) + 1
        local(i)  = size(value_scalar,i)
        stride(i) = 0
    end do

!   Populate test data based on process rank
    do j = lbound(value_scalar,1), ubound(value_scalar,1)
       value_scalar(j) = start(1) + j - 1
    end do
    do i = lbound(value_vector,1), ubound(value_vector,1)
        value_vector(i,:) = i*value_scalar(:)
    end do

!   Establish the parallel decomposition
    call esio_line_establish(h, global(1), start(1), local(1), ierr)
    ASSERT(ierr == 0)

!   Create a file
    call esio_file_create(h, filename, .false., ierr)
    ASSERT(ierr == 0)

!   Write a scalar-valued line using an explicitly-typed interface
    call esio_line_write_integer(h, "name_scalar", value_scalar,            &
                                 stride(1),                                 &
                                 "scalar-valued line", ierr)
    ASSERT(ierr == 0)

!   Overwrite the scalar-valued line using a generic interface
    call esio_line_write(h, "name_scalar", value_scalar,            &
                         stride(1),                                 &
                         "scalar-valued line", ierr)
    ASSERT(ierr == 0)

!   If contiguous, repeat the above using a named comment, no stride, no ierr
    if (stride(1) == 0) then
        call esio_line_write_integer(h, "name_scalar", value_scalar,        &
                                     comment="scalar-valued line")
        call esio_line_write        (h, "name_scalar", value_scalar,        &
                                     comment="scalar-valued line")
    end if

!   Write a vector-valued line using an explicitly-typed interface
    call esio_line_writev_integer(h, "name_vector", value_vector,           &
                                  ncomponents, stride(1),                   &
                                  "vector-valued line", ierr)
    ASSERT(ierr == 0)

!   Overwrite the vector-valued line using a generic interface
    call esio_line_writev(h, "name_vector", value_vector,   &
                          ncomponents, stride(1),           &
                          "vector-valued line", ierr)
    ASSERT(ierr == 0)

!   Close the file
    call esio_file_close(h, ierr)
    ASSERT(ierr == 0)

!   Re-open the file read-only
    call esio_file_open(h, filename, .false., ierr)
    ASSERT(ierr == 0)

!   Check the scalar-valued line's size and data
    call esio_line_size(h, "name_scalar", j, ierr)
    ASSERT(ierr == 0)
    ASSERT(j == global(1))
    call esio_line_size(h, "DOESNOTEXIST", j, ierr)
    ASSERT(ierr /= 0)
    ASSERT(j == global(1))  ! Preserved per feature #2419?
    call esio_line_sizev(h, "name_scalar", i, j, ierr)
    ASSERT(ierr == 0)
    ASSERT(i == 1)
    ASSERT(j == global(1))
    call esio_line_sizev(h, "DOESNOTEXIST", i, j, ierr)
    ASSERT(ierr /= 0)
    ASSERT(i == 1)          ! Preserved per feature #2419?
    ASSERT(j == global(1))  ! Preserved per feature #2419?
    call esio_line_read_integer(h, "name_scalar", buffer_scalar,          &
                                stride(1),                                &
                                ierr)
    ASSERT(ierr == 0)
    ASSERT(all(buffer_scalar == value_scalar))
    buffer_scalar = 0
    call esio_line_read(h, "name_scalar", buffer_scalar,          &
                        stride(1),                                &
                        ierr)
    ASSERT(ierr == 0)
    ASSERT(all(buffer_scalar == value_scalar))
    buffer_scalar = 0

!   Check the vector-valued line's size and data
    call esio_line_sizev(h, "name_vector", i, j, ierr)
    ASSERT(ierr == 0)
    ASSERT(i == ncomponents)
    ASSERT(j == global(1))
    call esio_line_readv_integer(h, "name_vector", buffer_vector,          &
                                 ncomponents,                              &
                                 stride(1),                                &
                                 ierr)
    ASSERT(ierr == 0)
    ASSERT(all(buffer_vector == value_vector))
    buffer_vector = 0
    call esio_line_readv(h, "name_vector", buffer_vector,          &
                         ncomponents,                              &
                         stride(1),                                &
                         ierr)
    ASSERT(ierr == 0)
    ASSERT(all(buffer_vector == value_vector))
    buffer_vector = 0

!   Close the file
    call esio_file_close(h, ierr)
    ASSERT(ierr == 0)

    call testframework_teardown()

end program line_int_f
