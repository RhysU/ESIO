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

program basic_f

    use mpi
    use esio
    use testframework

    implicit none

    logical :: file_exists

    call testframework_setup()

!   Check that Fortran thinks the right file does not (yet) exist
    inquire (file=trim(filename), exist=file_exists)
    ASSERT(.not. file_exists)

!   Open and close a unique file with overwrite disabled
    call esio_file_create(h, filename, .false., ierr)
    ASSERT(ierr == 0)
    call esio_file_close(h, ierr)
    ASSERT(ierr == 0)

!   Check that Fortran thinks the right file exists
    inquire (file=trim(filename), exist=file_exists)
    ASSERT(file_exists)

!   Open and close an old file with overwrite enabled
    call esio_file_create(h, filename, .true., ierr)
    ASSERT(ierr == 0)
    call esio_file_close(h, ierr)
    ASSERT(ierr == 0)

    call testframework_teardown()

end program basic_f
