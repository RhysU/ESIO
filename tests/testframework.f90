!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! esio 0.1.7: ExaScale IO library for turbulence simulation restart files
!! http://pecos.ices.utexas.edu/
!!
!! Copyright (C) 2010, 2011 The PECOS Development Team
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
  use ifport, only: abort
#endif /* __INTEL_COMPILER */
  use mpi
  use esio
  use esio_c_binding

  implicit none  ! Nothing implicit
  private        ! Everything default private

  integer,            public :: ierr, world_size, world_rank, output
  character(len=255), public :: input_dir, output_dir, filetemplate, filename
  integer,            public :: ndims, cart_comm, cart_rank
  type(esio_handle),  public :: h

  integer, public, allocatable, dimension(:) :: dims, coords
  integer, public, allocatable, dimension(:) :: global, start, local, stride
  logical, public                            :: verbose = .false.

  public :: testframework_setup, testframework_teardown
  public :: testframework_assert, testframework_unlink

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine testframework_setup (testsource, dimensionality)

    use, intrinsic :: iso_fortran_env, only: output_unit

    character(len=*), intent(in)           :: testsource
    integer,          intent(in), optional :: dimensionality
    logical, allocatable, dimension(:)     :: periods  ! OpenMPI workaround

!   Determine process topology dimensionality and allocate storage
    if (present(dimensionality)) then
      ndims = dimensionality
    else
      ndims = 0
    end if
    allocate ( dims(ndims), coords(ndims),                              &
               global(ndims), start(ndims), local(ndims), stride(ndims) )

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
    if (.not. create_testfiletemplate(output_dir, &
                                      testsource, &
                                      filetemplate)) then
        call MPI_Abort (MPI_COMM_WORLD, 1, ierr)
    end if
    if (world_rank == 0) then
      if (.not. create_testfilename(filetemplate, filename)) then
        call MPI_Abort (MPI_COMM_WORLD, 1, ierr)
      end if
      if (verbose) then
        write (output, *) "Generated test filename: ", trim(filename)
      end if
    end if
    call MPI_Bcast (filename, len(filename), MPI_CHARACTER,  &
                    0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) call MPI_Abort (MPI_COMM_WORLD, 1, ierr)

!   Create an ndims-dimensional process topology on cart_comm when ndims > 0
!   Otherwise duplicate MPI_COMM_WORLD details into cart_comm, cart_rank
    if (ndims > 0) then
        dims(:) = 0
        call MPI_Dims_create (world_size, ndims, dims, ierr)
        if (verbose) then
          write (output, *) "Test topology is [", dims, "] on ", &
                            world_size, " ranks"
        end if
        if (ierr /= MPI_SUCCESS) call MPI_Abort (MPI_COMM_WORLD, 1, ierr)
        allocate ( periods(ndims) )
        periods(:) = .false.
        call MPI_Cart_create (MPI_COMM_WORLD, ndims, dims, periods, .true., &
                              cart_comm, ierr)
        deallocate ( periods )
        if (ierr /= MPI_SUCCESS) call MPI_Abort (MPI_COMM_WORLD, 1, ierr)
    else
        call MPI_Comm_dup (MPI_COMM_WORLD, cart_comm, ierr)
        if (ierr /= MPI_SUCCESS) call MPI_Abort (MPI_COMM_WORLD, 1, ierr)
    end if
    CALL MPI_Comm_rank (cart_comm, cart_rank, ierr)
    if (ierr /= MPI_SUCCESS) call MPI_Abort (MPI_COMM_WORLD, 1, ierr)
    if (ndims > 0) then
        CALL MPI_Cart_coords (cart_comm, cart_rank, ndims, coords, ierr)
        if (ierr /= MPI_SUCCESS) call MPI_Abort (MPI_COMM_WORLD, 1, ierr)
    end if

!   Initialize an ESIO handle against the new communicator
    call esio_handle_initialize (h, cart_comm, ierr)
    if (ierr /= MPI_SUCCESS) call MPI_Abort (MPI_COMM_WORLD, 1, ierr)

  end subroutine testframework_setup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine testframework_teardown ()

!   Finalize the ESIO handle
    call esio_handle_finalize (h, ierr)
    if (ierr /= MPI_SUCCESS) call MPI_Abort (MPI_COMM_WORLD, 1, ierr)

!   Finalize the communicator
    call MPI_Comm_free (cart_comm, ierr)
    if (ierr /= MPI_SUCCESS) call MPI_Abort (MPI_COMM_WORLD, 1, ierr)

!   Deallocate working storage for process topology
    deallocate ( dims, coords, global, start, local, stride )

!   Close rank-dependent output unit
    if (world_rank /= 0) then
      close (7)
    end if

!   Attempt to delete the named temporary file, if it exists
    call testframework_unlink(filename)

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

! Really?!  http://software.intel.com/en-us/forums/showthread.php?t=66887
#ifndef __INTEL_COMPILER
    flush (output_unit)
    flush (error_unit)
#endif
    if (verbose .and. condition) then
      write (output_unit, *) "Assertion passed at ", sourcefile, ":", line, &
                            " on rank ", world_rank
    else if (.not. condition) then
      write (error_unit, *) "Assertion failed at ", sourcefile, ":", line,  &
                            " on rank ", world_rank
      call MPI_Abort (MPI_COMM_WORLD, 1, ierr)
      stop
    end if
#ifndef __INTEL_COMPILER
    flush (error_unit)
    flush (output_unit)
#endif

  end subroutine testframework_assert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function create_testfiletemplate(dir, filename, filetemplate)  &
           result (success)

    use, intrinsic :: iso_c_binding, only: c_ptr, c_char

    logical                         :: success
    character(len=*), intent(in)    :: dir
    character(len=*), intent(in)    :: filename
    character(len=*), intent(inout) :: filetemplate
    type(c_ptr)                     :: tmp_p

    interface
      function impl (dir, filename) bind (C, name="create_testfiletemplate")
        import
        type(c_ptr)                              :: impl
        character(len=1,kind=c_char), intent(in) :: dir(*)
        character(len=1,kind=c_char), intent(in) :: filename(*)
      end function impl
    end interface

    success = .false.
    tmp_p = impl(esio_f_c_string(dir), esio_f_c_string(filename))
    if (esio_c_f_stringcopy(tmp_p, filetemplate)) then
      success = .true.
    end if
    call esio_c_free(tmp_p)

  end function create_testfiletemplate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function create_testfilename(filetemplate, filename) result (success)

    use, intrinsic :: iso_c_binding, only: c_ptr, c_char

    logical                         :: success
    character(len=*), intent(in)    :: filetemplate
    character(len=*), intent(inout) :: filename
    type(c_ptr)                     :: tmp_p

    interface
      function impl (filetemplate) bind (C, name="create_testfilename")
        import
        type(c_ptr)                              :: impl
        character(len=1,kind=c_char), intent(in) :: filetemplate(*)
      end function impl
    end interface

    success = .false.
    tmp_p = impl(esio_f_c_string(filetemplate))
    if (esio_c_f_stringcopy(tmp_p, filename)) then
      success = .true.
    end if
    call esio_c_free(tmp_p)

  end function create_testfilename

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine testframework_unlink (filename)

#ifdef __INTEL_COMPILER
  use ifport, only: unlink
#endif /* __INTEL_COMPILER */

    character(len=*), intent(in) :: filename
    logical                      :: file_exists

!   Attempt to delete the named temporary file, if it exists
    if (.not. verbose .and. world_rank == 0) then
      inquire (file=trim(filename), exist=file_exists)
      if (file_exists) ierr = unlink(trim(filename))
    end if

  end subroutine testframework_unlink

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module testframework
