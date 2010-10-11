!!--------------------------------------------------------------------------
!!--------------------------------------------------------------------------
!!
!! ESIO - an IO library for scientific computing
!!
!! Copyright (C) 2010 The PECOS Development Team
!!
!! ESIO is free software; you can redistribute it and/or
!! modify it under the terms of the Version 2.1 GNU Lesser General
!! Public License as published by the Free Software Foundation.
!!
!! ESIO is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!! Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public
!! License along with ESIO; if not, write to the Free Software
!! Foundation, Inc. 51 Franklin Street, Fifth Floor,
!! Boston, MA  02110-1301  USA
!!
!!--------------------------------------------------------------------------
!!
!! esio.f90: Fortran module describing ESIO's public Fortran API
!!
!! $Id$
!!--------------------------------------------------------------------------
!!--------------------------------------------------------------------------

! TODO Document Fortran API

module esio

! Use ISO_C_BINDING for interoperation with ESIO's C implementation.
! Rename C_PTR to ESIO_STATE to make the handle object types more opaque.
! State traverses the language boundary as TYPE(ESIO_STATE)s.
! Renaming also prevents collision if caller uses ISO_C_BINDING directly.
! Other interoperation details are kept hidden from the client.
  use, intrinsic :: iso_c_binding, only: c_char,             &
                                         c_double,           &
                                         c_double_complex,   &
                                         c_float,            &
                                         c_float_complex,    &
                                         c_int,              &
                                         c_null_char,        &
                                         esio_state => c_ptr
  private :: c_char, c_double, c_double_complex
  private :: c_float, c_float_complex, c_int, c_null_char
  public  :: esio_state

! Public functionality
  public :: esio_init, esio_finalize
  public :: esio_file_create, esio_file_open, esio_file_close
  public :: esio_field_write, esio_field_write_double, esio_field_write_single

! Generic, precision-agnostic field interfaces
  interface esio_field_write
    module procedure esio_field_write_double, esio_field_write_single
  end interface

! Private utilities to convert types from Fortran- to C-style.
  private :: f_c_string
  private :: f_c_logical

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(esio_state) function esio_init (comm)

    integer, intent(in) :: comm

    ! See C routine esio_init_fortran re: MPI communicator interoperation
    interface
      import
      function impl (comm) bind (C, name="esio_init_fortran")
        type(esio_state)           :: impl
        integer, intent(in), value :: comm  ! Note integer not integer(c_int)
      end function impl
    end interface

    esio_init = impl(comm)

  end function esio_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_finalize ()

    interface
      import
      function impl () bind (C, name="esio_finalize")
        integer(c_int) :: impl
      end function impl
    end interface

    esio_finalize = impl()

  end function esio_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_file_create (state, file, overwrite)

    type(esio_state), intent(in) :: state
    character(len=*), intent(in) :: file
    logical,          intent(in) :: overwrite

    interface
      import
      function impl (state, file, overwrite) bind (C, name="esio_file_create")
        integer(c_int)                                  :: impl
        type(esio_state),             intent(in), value :: state
        character(len=1,kind=c_char), intent(in), value :: file
        logical,                      intent(in), value :: overwrite
      end function impl
    end interface

    esio_file_create = impl(state, f_c_string(file), f_c_logical(overwrite))

  end function esio_file_create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_file_open (state, file, readwrite)

    type(esio_state), intent(in) :: state
    character(len=*), intent(in) :: file
    logical,          intent(in) :: readwrite

    interface
      import
      function impl (state, file, readwrite) bind (C, name="esio_file_open")
        integer(c_int)                                  :: impl
        type(esio_state),             intent(in), value :: state
        character(len=1,kind=c_char), intent(in), value :: file
        logical,                      intent(in), value :: readwrite
      end function impl
    end interface

    esio_file_open = impl(state, f_c_string(file), f_c_logical(readwrite))

  end function esio_file_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_file_close (state)

    type(esio_state), intent(in) :: state

    interface
      import
      function impl (state) bind (C, name="esio_file_close")
        integer(c_int)                      :: impl
        type(esio_state), intent(in), value :: state
      end function impl
    end interface

    esio_file_close = impl(state, f_c_string(file), f_c_logical(readwrite))

  end function esio_file_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! TODO Use Autoconf to find c_double value at compile time and use in API.
! Would eliminate "leaking" c_double in the API and leave, e.g., real(8).

  integer function esio_field_write_double (state, name, data, &
                                            na, ast, asz,      &
                                            nb, bst, bsz,      &
                                            nc, cst, csz)

    type(esio_state), intent(in) :: state
    character(len=*), intent(in) :: name
    real(c_double),   intent(in) :: data(*)
    integer,          intent(in) :: na, ast, asz
    integer,          intent(in) :: nb, bst, bsz
    integer,          intent(in) :: nc, cst, csz

    interface
      import
      function impl (state, name, data, &
                     na, ast, asz,      &
                     nb, bst, bsz,      &
                     nc, cst, csz) bind (C, name="esio_field_write_double")
        integer(c_int)                                  :: impl
        type(esio_state),             intent(in), value :: state
        character(len=1,kind=c_char), intent(in), value :: name
        real(c_double),               intent(in), value :: data(*)
        integer(c_int),               intent(in), value :: na, ast, asz
        integer(c_int),               intent(in), value :: nb, bst, bsz
        integer(c_int),               intent(in), value :: nc, cst, csz
      end function impl
    end interface

!   Note conversion from one- to zero-based starting offsets
    esio_field_write_double = impl(state, f_c_string(name), data, &
                                   na, ast - 1, asz,              &
                                   nb, bst - 1, bsz,              &
                                   nc, cst - 1, csz)

  end function esio_field_write_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! TODO Use Autoconf to find c_float value at compile time and use in API.
! Would eliminate "leaking" c_float in the API and leave, e.g., real(4).

  integer function esio_field_write_single (state, name, data, &
                                            na, ast, asz,      &
                                            nb, bst, bsz,      &
                                            nc, cst, csz)

    type(esio_state), intent(in) :: state
    character(len=*), intent(in) :: name
    real(c_float),    intent(in) :: data(*)
    integer,          intent(in) :: na, ast, asz
    integer,          intent(in) :: nb, bst, bsz
    integer,          intent(in) :: nc, cst, csz

    interface
      import
      function impl (state, name, data, &
                     na, ast, asz,      &
                     nb, bst, bsz,      &
                     nc, cst, csz) bind (C, name="esio_field_write_float")
        integer(c_int)                                  :: impl
        type(esio_state),             intent(in), value :: state
        character(len=1,kind=c_char), intent(in), value :: name
        real(c_float),                intent(in), value :: data(*)
        integer(c_int),               intent(in), value :: na, ast, asz
        integer(c_int),               intent(in), value :: nb, bst, bsz
        integer(c_int),               intent(in), value :: nc, cst, csz
      end function impl
    end interface

!   Note conversion from one- to zero-based starting offsets
    esio_field_write_single = impl(state, f_c_string(name), data, &
                                   na, ast - 1, asz,              &
                                   nb, bst - 1, bsz,              &
                                   nc, cst - 1, csz)

  end function esio_field_write_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Convert Fortran-style string to a C-style one
! See http://fortranwiki.org/fortran/show/Generating+C+Interfaces#strings
  pure function f_c_string (f_string) result (c_string)

    character(len=*), intent(in) :: f_string
    character(len=1,kind=c_char) :: c_string(len_trim(f_string)+1)

    c_string = trim(f_string) // c_null_char

  end function f_c_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Convert Fortran LOGICAL to a C-style true/false INT
  elemental function f_c_logical (f_logical) result (c_logical)

    logical, intent(in) :: f_logical
    integer(c_int)      :: c_logical

    if (f_logical) then
      c_logical = 1_c_int
    else
      c_logical = 0_c_int
    end if

  end function f_c_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module esio
