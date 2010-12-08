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

program basic_f

    use mpi
    use esio
    use testframework

    implicit none

    logical :: file_exists

    call testframework_setup()

!   Check that Fortran thinks the right file does not (yet) exist
    inquire (file=trim(filename), exist=file_exists)
    if (file_exists) call MPI_Abort(MPI_COMM_WORLD, 1, ierr)

!   Open and close a unique file with overwrite disabled
    call esio_file_create(h, filename, .false., ierr)
    if (ierr /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    call esio_file_close(h, ierr)
    if (ierr /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD, 1, ierr)

!   Check that Fortran thinks the right file exists
    inquire (file=trim(filename), exist=file_exists)
    if (.not. file_exists) call MPI_Abort(MPI_COMM_WORLD, 1, ierr)

!   Open and close an old file with overwrite enabled
    call esio_file_create(h, filename, .true., ierr)
    if (ierr /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    call esio_file_close(h, ierr)
    if (ierr /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD, 1, ierr)

    call testframework_teardown()

end program basic_f
