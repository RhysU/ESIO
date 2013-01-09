!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! ESIO 0.1.8: ExaScale IO library for turbulence simulation restart files
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

program basic_f

    use esio
    use testframework

    implicit none

    integer            :: size, rank
    logical            :: file_exists
    character(len=256) :: file_path, template, restart0, restart1

!   Used only for checking establish/established
    integer, parameter :: cg = 16, cs = 7, cl = 6
    integer, parameter :: bg = 13, bs = 4, bl = 5
    integer, parameter :: ag = 11, as = 2, al = 3
    integer            :: tmp_cg = 555, tmp_cs, tmp_cl
    integer            :: tmp_bg = 555, tmp_bs, tmp_bl
    integer            :: tmp_ag = 555, tmp_as, tmp_al

!   Assert status codes are publicly exposed in the API
!   Refer to feature #2184 and bug #2410 for why ifdefs appear.
#ifndef __INTEL_COMPILER
#if defined(__GNUC__) && (__GNUC__ > 4 || __GNUC_MINOR__ > 5)
    ASSERT(ESIO_SUCCESS  ==  0)
    ASSERT(ESIO_EFAULT   ==  3)
    ASSERT(ESIO_EINVAL   ==  4)
    ASSERT(ESIO_EFAILED  ==  5)
    ASSERT(ESIO_ESANITY  ==  7)
    ASSERT(ESIO_ENOMEM   ==  8)
    ASSERT(ESIO_NOTFOUND ==  9)
#endif
#endif

    call testframework_setup(__FILE__)

!   Check that the handle correctly reports world size and world rank
    call esio_handle_comm_size(h, size, ierr)
    ASSERT(ierr == 0)
    ASSERT(size == world_size)
    call esio_handle_comm_rank(h, rank, ierr)
    ASSERT(ierr == 0)
    ASSERT(rank == world_rank)

!   Check that we can establish and retrieve line decomposition details
    call esio_line_established(h, tmp_ag, tmp_as, tmp_al, ierr)
    ASSERT(ierr == 0)
    ASSERT(tmp_ag == 0)
    ASSERT(tmp_as == 1)  ! One-based
    ASSERT(tmp_al == 0)
    call esio_line_establish(h, ag, as, al, ierr)
    ASSERT(ierr == 0)
    call esio_line_established(h, tmp_ag, tmp_as, tmp_al, ierr)
    ASSERT(ierr == 0)
    ASSERT(tmp_ag == ag)
    ASSERT(tmp_as == as)
    ASSERT(tmp_al == al)

!   Check that we can establish and retrieve plane decomposition details
    call esio_plane_established(h, tmp_ag, tmp_as, tmp_al,       &
                                   tmp_bg, tmp_bs, tmp_bl, ierr)
    ASSERT(ierr == 0)
    ASSERT(tmp_ag == 0)
    ASSERT(tmp_as == 1)  ! One-based
    ASSERT(tmp_al == 0)
    ASSERT(tmp_bg == 0)
    ASSERT(tmp_bs == 1)  ! One-based
    ASSERT(tmp_bl == 0)
    call esio_plane_establish(h, ag, as, al,       &
                                 bg, bs, bl, ierr)
    ASSERT(ierr == 0)
    call esio_plane_established(h, tmp_ag, tmp_as, tmp_al,       &
                                   tmp_bg, tmp_bs, tmp_bl, ierr)
    ASSERT(ierr == 0)
    ASSERT(tmp_ag == ag)
    ASSERT(tmp_as == as)
    ASSERT(tmp_al == al)
    ASSERT(tmp_bg == bg)
    ASSERT(tmp_bs == bs)
    ASSERT(tmp_bl == bl)

!   Check that we can establish and retrieve field decomposition details
    call esio_field_established(h, tmp_ag, tmp_as, tmp_al,       &
                                   tmp_bg, tmp_bs, tmp_bl,       &
                                   tmp_cg, tmp_cs, tmp_cl, ierr)
    ASSERT(ierr == 0)
    ASSERT(tmp_ag == 0)
    ASSERT(tmp_as == 1)  ! One-based
    ASSERT(tmp_al == 0)
    ASSERT(tmp_bg == 0)
    ASSERT(tmp_bs == 1)  ! One-based
    ASSERT(tmp_bl == 0)
    ASSERT(tmp_cg == 0)
    ASSERT(tmp_cs == 1)  ! One-based
    ASSERT(tmp_cl == 0)
    call esio_field_establish(h, ag, as, al,       &
                                 bg, bs, bl,       &
                                 cg, cs, cl, ierr)
    ASSERT(ierr == 0)
    call esio_field_established(h, tmp_ag, tmp_as, tmp_al,       &
                                   tmp_bg, tmp_bs, tmp_bl,       &
                                   tmp_cg, tmp_cs, tmp_cl, ierr)
    ASSERT(ierr == 0)
    ASSERT(tmp_ag == ag)
    ASSERT(tmp_as == as)
    ASSERT(tmp_al == al)
    ASSERT(tmp_bg == bg)
    ASSERT(tmp_bs == bs)
    ASSERT(tmp_bl == bl)
    ASSERT(tmp_cg == cg)
    ASSERT(tmp_cs == cs)
    ASSERT(tmp_cl == cl)

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
    call testframework_unlink(restart0)
    call testframework_unlink(restart1)

    call testframework_teardown()

end program basic_f
