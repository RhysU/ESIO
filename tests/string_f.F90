!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! ESIO 0.2.0: ExaScale IO library for turbulence simulation restart files
!! http://red.ices.utexas.edu/projects/esio/
!!
!! Copyright (C) 2010-2014 The PECOS Development Team
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

program string_f

    use esio
    use testframework

    implicit none

    character(len=*), parameter :: value = "testTESTDATAdata"
    character(len=255)          :: buffer

    call testframework_setup(__FILE__)

!   Create a file, write a string, and close it
    call esio_file_create(h, filename, .false., ierr)
    ASSERT(ierr == 0)
    call esio_string_set(h, "/", "name", value, ierr)
    ASSERT(ierr == 0)
    call esio_file_close(h, ierr)
    ASSERT(ierr == 0)

!   Re-open the file read-only, check the string, and close it
    call esio_file_open(h, filename, .false., ierr)
    ASSERT(ierr == 0)
    call esio_string_get(h, "/", "name", buffer, ierr)
    ASSERT(ierr == 0)
    ASSERT(value == trim(buffer))
    call esio_file_close(h, ierr)
    ASSERT(ierr == 0)

    call testframework_teardown()

end program string_f
