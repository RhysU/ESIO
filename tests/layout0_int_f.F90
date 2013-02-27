!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! ESIO 0.1.9: ExaScale IO library for turbulence simulation restart files
!! http://red.ices.utexas.edu/projects/esio/
!!
!! Copyright (C) 2010, 2011, 2012, 2013 The PECOS Development Team
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

program layout0_int_f

    use esio
    use testframework

    implicit none

    integer            :: value_scalar(3,4,5),   buffer_scalar(3,4,5)
    integer            :: value_vector(2,3,4,5), buffer_vector(2,3,4,5)
    integer, parameter :: ncomponents = size(value_vector, 1)
    integer            :: i, j, k, l

    ! Three-dimensional Cartesian topology
    call testframework_setup (__FILE__, 3)

!   Set and retrieve the field layout to sniff test the API
    call esio_field_layout_count(i, ierr)
    ASSERT(ierr == 0)
    ASSERT(i > 0)
    call esio_field_layout_set(h, 1, ierr)
    ASSERT(ierr == 0)
    call esio_field_layout_get(h, i, ierr)
    ASSERT(ierr == 0)
    ASSERT(i == 1)

!   Set layout 0
    call esio_field_layout_set(h, 0, ierr)
    ASSERT(ierr == 0)
    call esio_field_layout_get(h, i, ierr)
    ASSERT(ierr == 0)
    ASSERT(i == 0)

!   Prepare homogeneous test data distribution across the topology
    do i = 1, ndims
        global(i) = size(value_scalar,i) * dims(i)
        start(i)  = size(value_scalar,i) * coords(i) + 1
        local(i)  = size(value_scalar,i)
        stride(i) = 0
    end do

!   Populate test data based on process rank
    do l = lbound(value_scalar,3), ubound(value_scalar,3)
        do k = lbound(value_scalar,2), ubound(value_scalar,2)
            do j = lbound(value_scalar,1), ubound(value_scalar,1)
                value_scalar(j,k,l) = (start(1) + j - 1) &
                                    * (start(2) + k - 1) &
                                    * (start(3) + l - 1)
            end do
        end do
    end do
    do i = lbound(value_vector,1), ubound(value_vector,1)
        value_vector(i,:,:,:) = i*value_scalar(:,:,:)
    end do

!   Establish the parallel decomposition
    call esio_field_establish(h, global(1), start(1), local(1),       &
                                 global(2), start(2), local(2),       &
                                 global(3), start(3), local(3), ierr)

!   Create a file
    call esio_file_create(h, filename, .false., ierr)
    ASSERT(ierr == 0)

!   Write a scalar-valued field using an explicitly-typed interface
    call esio_field_write_integer(h, "name_scalar", value_scalar,           &
                                  stride(1), stride(2), stride(3),          &
                                  "scalar-valued field", ierr)
    ASSERT(ierr == 0)

!   Overwrite the scalar-valued field using a generic interface
    call esio_field_write(h, "name_scalar", value_scalar,           &
                          stride(1), stride(2), stride(3),          &
                          "scalar-valued field", ierr)
    ASSERT(ierr == 0)

!   Write a vector-valued field using an explicitly-typed interface
    call esio_field_writev_integer(h, "name_vector", value_vector,           &
                                   ncomponents,                              &
                                   stride(1), stride(2), stride(3),          &
                                   "vector-valued field", ierr)
    ASSERT(ierr == 0)

!   Overwrite the vector-valued field using a generic interface
    call esio_field_writev(h, "name_vector", value_vector,           &
                           ncomponents,                              &
                           stride(1), stride(2), stride(3),          &
                           "vector-valued field", ierr)
    ASSERT(ierr == 0)

!   If contiguous, repeat the above using a named comment, no stride, no ierr
    if (stride(1) == 0 .and. stride(2) == 0 .and. stride(3) == 0) then
        call esio_field_writev_integer(h, "name_vector", value_vector, &
                                       ncomponents,                    &
                                       comment="vector-valued field")
        call esio_field_writev        (h, "name_vector", value_vector, &
                                       ncomponents,                    &
                                       comment="vector-valued field")
    end if

!   Close the file
    call esio_file_close(h, ierr)
    ASSERT(ierr == 0)

!   Re-open the file read-only
    call esio_file_open(h, filename, .false., ierr)
    ASSERT(ierr == 0)

!   Check the scalar-valued field's size and data
    call esio_field_size(h, "name_scalar", j, k, l, ierr)
    ASSERT(ierr == 0)
    ASSERT(j == global(1))
    ASSERT(k == global(2))
    ASSERT(l == global(3))
    call esio_field_size(h, "DOESNOTEXIST", j, k, l, ierr)
    ASSERT(ierr /= 0)
    ASSERT(j == global(1))  ! Preserved per feature #2419?
    ASSERT(k == global(2))  ! Preserved per feature #2419?
    ASSERT(l == global(3))  ! Preserved per feature #2419?
    call esio_field_sizev(h, "name_scalar", i, j, k, l, ierr)
    ASSERT(ierr == 0)
    ASSERT(i == 1)
    ASSERT(j == global(1))
    ASSERT(k == global(2))
    ASSERT(l == global(3))
    call esio_field_sizev(h, "DOESNOTEXIST", i, j, k, l, ierr)
    ASSERT(ierr /= 0)
    ASSERT(i == 1)          ! Preserved per feature #2419?
    ASSERT(j == global(1))  ! Preserved per feature #2419?
    ASSERT(k == global(2))  ! Preserved per feature #2419?
    ASSERT(l == global(3))  ! Preserved per feature #2419?
    call esio_field_read_integer(h, "name_scalar", buffer_scalar,          &
                                 stride(1), stride(2), stride(3),          &
                                 ierr)
    ASSERT(ierr == 0)
    ASSERT(all(buffer_scalar == value_scalar))
    buffer_scalar = 0
    call esio_field_read(h, "name_scalar", buffer_scalar,          &
                         stride(1), stride(2), stride(3),          &
                         ierr)
    ASSERT(ierr == 0)
    ASSERT(all(buffer_scalar == value_scalar))
    buffer_scalar = 0

!   Check the vector-valued field's size and data
    call esio_field_sizev(h, "name_vector", i, j, k, l, ierr)
    ASSERT(ierr == 0)
    ASSERT(i == ncomponents)
    ASSERT(j == global(1))
    ASSERT(k == global(2))
    ASSERT(l == global(3))
    call esio_field_readv_integer(h, "name_vector", buffer_vector,          &
                                  ncomponents,                              &
                                  stride(1), stride(2), stride(3),          &
                                  ierr)
    ASSERT(ierr == 0)
    ASSERT(all(buffer_vector == value_vector))
    buffer_vector = 0
    call esio_field_readv(h, "name_vector", buffer_vector,          &
                          ncomponents,                              &
                          stride(1), stride(2), stride(3),          &
                          ierr)
    ASSERT(ierr == 0)
    ASSERT(all(buffer_vector == value_vector))
    buffer_vector = 0

!   Close the file
    call esio_file_close(h, ierr)
    ASSERT(ierr == 0)

    call testframework_teardown()

end program layout0_int_f
