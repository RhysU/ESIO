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

program string_f

    use mpi
    use esio
    use testframework

    implicit none

    character(len=*), parameter :: value = "testTESTDATAdata"
    character(len=255)          :: buffer

    call testframework_setup()

!   Create a file, write a string, and close it
    call esio_file_create(h, filename, .false., ierr)
    if (ierr /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    call esio_string_set(h, "name", value, ierr)
    if (ierr /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    call esio_file_close(h, ierr)
    if (ierr /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD, 1, ierr)

!   Re-open the file read-only, check the string, and close it
    call esio_file_open(h, filename, .false., ierr)
    if (ierr /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    call esio_string_get(h, "name", buffer, ierr)
    if (ierr /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    if (value /= trim(buffer)) call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    call esio_file_close(h, ierr)
    if (ierr /= MPI_SUCCESS) call MPI_Abort(MPI_COMM_WORLD, 1, ierr)

    call testframework_teardown()

end program string_f
