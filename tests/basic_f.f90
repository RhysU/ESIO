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

    use esio
    use testframework

    implicit none

    logical            :: file_exists
    character(len=256) :: file_path, template, restart0, restart1

    call testframework_setup(__FILE__)

!   Check that Fortran thinks the right file does not (yet) exist
    inquire (file=trim(filename), exist=file_exists)
    ASSERT(.not. file_exists)

!   Open and close a unique file with overwrite disabled
    call esio_file_create(h, filename, .false., ierr)
    ASSERT(ierr == 0)
    call esio_file_path(h, file_path, ierr)
    ASSERT(ierr == 0)
    call esio_file_close(h, ierr)
    ASSERT(ierr == 0)

!   Check that Fortran thinks the right file exists
    inquire (file=trim(filename), exist=file_exists)
    ASSERT(file_exists)

!   Check that Fortran thinks the absolute file path exists
    inquire (file=trim(file_path), exist=file_exists)
    ASSERT(file_exists)
!   Check that the absolute file path begins with a '/'
    ASSERT(scan(file_path, "/") == 1)

!   Open, flush, and close an old file with overwrite enabled
    call esio_file_create(h, filename, .true., ierr)
    ASSERT(ierr == 0)
    call esio_file_flush(h, ierr)
    ASSERT(ierr == 0)
    call esio_file_close(h, ierr)
    ASSERT(ierr == 0)

!   Prepare to check the restart convenience methods
    template = trim(filename) // "#"
    restart0 = trim(filename) // "0"
    restart1 = trim(filename) // "1"

!   Create with overwrite followed by esio_file_close_restart
    call esio_file_create(h, filename, .true., ierr)
    ASSERT(ierr == 0)
    call esio_file_close_restart(h, template, 2, ierr)
    ASSERT(ierr == 0)
    inquire (file=trim(filename), exist=file_exists)
    ASSERT(.not. file_exists)
    inquire (file=trim(restart0), exist=file_exists)
    ASSERT(      file_exists)
    inquire (file=trim(restart1), exist=file_exists)
    ASSERT(.not. file_exists)

!   Create without overwrite followed by esio_file_close_restart
    call esio_file_create(h, filename, .false., ierr)
    ASSERT(ierr == 0)
    call esio_file_close_restart(h, template, 2, ierr)
    ASSERT(ierr == 0)
    inquire (file=trim(filename), exist=file_exists)
    ASSERT(.not. file_exists)
    inquire (file=trim(restart0), exist=file_exists)
    ASSERT(      file_exists)
    inquire (file=trim(restart1), exist=file_exists)
    ASSERT(      file_exists)

!   Clean up restart files leftover from test
    if (.not. verbose .and. world_rank == 0) then
      inquire (file=trim(restart0), exist=file_exists)
      if (file_exists) ierr = unlink(trim(restart0))
      inquire (file=trim(restart1), exist=file_exists)
      if (file_exists) ierr = unlink(trim(restart1))
    end if

    call testframework_teardown()

end program basic_f
