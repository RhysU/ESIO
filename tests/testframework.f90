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

module testframework

#ifdef __INTEL_COMPILER
  use ifport, only: abort, unlink
#endif /* __INTEL_COMPILER */
  use mpi
  use esio
  use esio_c_binding

  implicit none  ! Nothing implicit
  private        ! Everything default private

  integer,            public :: ierr, world_size, world_rank, output
  character(len=255), public :: input_dir, output_dir, filename
  type(esio_handle),  public :: h

  logical :: verbose = .false.

  public :: testframework_setup, testframework_teardown, testframework_assert

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine testframework_setup ()

    use, intrinsic :: iso_fortran_env, only: output_unit

!   Initialize MPI
    call MPI_Init (ierr)
    if (ierr /= MPI_SUCCESS) call abort()
    call MPI_Comm_size (MPI_COMM_WORLD, world_size, ierr)
    if (ierr /= MPI_SUCCESS) call MPI_Abort (MPI_COMM_WORLD, 1, ierr)
    call MPI_Comm_rank (MPI_COMM_WORLD, world_rank, ierr)
    if (ierr /= MPI_SUCCESS) call MPI_Abort (MPI_COMM_WORLD, 1, ierr)

!   Initialize a rank-dependent output unit for progress messages
    if (world_rank == 0) then
      output = output_unit
    else
      output = 7
      open (7, file = '/dev/null', action = 'write')
    end if

!   Check if assertions should be verbose in nature
!   Uses input_dir as a scratchpad
    input_dir = ''
    call get_environment_variable("ESIO_TEST_VERBOSE", input_dir)
    verbose = len_trim(input_dir) > 0

!   Initialize a test-specific temporary filename
    call get_environment_variable("ESIO_TEST_INPUT_DIR",  input_dir)
    call get_environment_variable("ESIO_TEST_OUTPUT_DIR", output_dir)
    if (world_rank == 0) then
      if (.not. f_tempnam(output_dir, "etst", filename)) then
        call MPI_Abort (MPI_COMM_WORLD, 1, ierr)
      end if
      if (verbose) then
        write (output, *) "Generated test filename: ", trim(filename)
      end if
    end if
    call MPI_Bcast (filename, len(filename), MPI_CHARACTER,  &
                    0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) call MPI_Abort (MPI_COMM_WORLD, 1, ierr)

!   Initialize an ESIO handle against MPI_COMM_WORLD
    call esio_handle_initialize (h, MPI_COMM_WORLD)

  end subroutine testframework_setup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine testframework_teardown ()

    logical :: file_exists

!   Finalize the ESIO handle
    call esio_handle_finalize (h)

!   Close rank-dependent output unit
    if (world_rank /= 0) then
      close (7)
    end if

!   Attempt to delete the named temporary file, if it exists
    if (.not. verbose .and. world_rank == 0) then
      inquire (file=trim(filename), exist=file_exists)
      if (file_exists) ierr = unlink(trim(filename))
    end if

!   Finalize MPI
    call MPI_Finalize (ierr)
    if (ierr /= MPI_SUCCESS) call abort()

  end subroutine testframework_teardown

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine testframework_assert (condition, sourcefile, line)

    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit

    logical,          intent(in) :: condition
    character(len=*), intent(in) :: sourcefile
    integer,          intent(in) :: line

    flush (output_unit)
    flush (error_unit)
    if (verbose .and. condition) then
      write (output_unit, *) "Assertion passed at ", sourcefile, ":", line, &
                            " on rank ", world_rank
    else if (.not. condition) then
      write (error_unit, *) "Assertion failed at ", sourcefile, ":", line,  &
                            " on rank ", world_rank
      call MPI_Abort (MPI_COMM_WORLD, 1, ierr)
      stop
    end if
    flush (error_unit)
    flush (output_unit)

  end subroutine testframework_assert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function f_tempnam (dir, pfx, tempnam) result (success)

    use, intrinsic :: iso_c_binding, only: c_ptr, c_char

    logical                         :: success
    character(len=*), intent(in)    :: dir
    character(len=*), intent(in)    :: pfx
    character(len=*), intent(inout) :: tempnam
    type(c_ptr)                     :: tmp_p

!   See 'man 3 tempnam' for details on the C API
    interface
      function impl (dir, pfx) bind (C, name="tempnam")
        import
        type(c_ptr)                              :: impl
        character(len=1,kind=c_char), intent(in) :: dir(*)
        character(len=1,kind=c_char), intent(in) :: pfx(*)
      end function impl
    end interface

    success = .false.
    tmp_p = impl(esio_f_c_string(dir), esio_f_c_string(pfx))
    if (esio_c_f_stringcopy(tmp_p, tempnam)) then
      success = .true.
    end if
    call esio_c_free(tmp_p)

  end function f_tempnam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module testframework
