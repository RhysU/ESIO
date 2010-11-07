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

! TODO Document Fortran API
! TODO Allow Fortran to detect invalid handle before other failure
! TODO Allow Fortran to use customizable error handling
! TODO Disable error handling when ierr is present??

!> \file
!! Provides ESIO's Fortran-based public API following the library's
!! \ref concepts "usage concepts".  All methods in this header
!! invoke ESIO's \ref error.h "error handling mechanisms" on failure.

module esio

! Use ISO_C_BINDING for interoperation with ESIO's C implementation.
! Rename C_PTR to ESIO_HANDLE to make the handle object types more opaque.
! State traverses the language boundary as TYPE(ESIO_HANDLE)s.
! Renaming also prevents collision if caller uses ISO_C_BINDING directly.
  use, intrinsic :: iso_c_binding, only: c_char,             &
                                         c_double,           &
                                         c_double_complex,   &
                                         c_float,            &
                                         c_float_complex,    &
                                         c_int,              &
                                         c_null_char,        &
                                         c_associated,       &
                                         c_funptr,           &
                                         esio_handle => c_ptr

  implicit none  ! Nothing implicit

! C interoperation details are kept hidden from the client...
  private :: c_char, c_double, c_double_complex
  private :: c_float, c_float_complex, c_int, c_null_char
  private :: c_associated, c_funptr
! ... with the exception of our opaque handle object type.
  public :: esio_handle

! Generic, precision-agnostic interfaces atop the public API
  interface esio_field_write
    module procedure esio_field_write_double
    module procedure esio_field_write_single
  end interface

  interface esio_field_read
    module procedure esio_field_read_double
    module procedure esio_field_read_single
  end interface

! Internal helper routines
  private :: f_c_string, f_c_logical

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_initialize (h, comm, ierr)

    type(esio_handle), intent(out)           :: h
    integer,           intent(in)            :: comm
    integer,           intent(out), optional :: ierr

    ! See C routine esio_initialize_fortran re: MPI communicator interoperation
    interface
      function impl (comm) bind (C, name="esio_initialize_fortran")
        import
        type(esio_handle)          :: impl
        integer, intent(in), value :: comm  ! Note integer not integer(c_int)
      end function impl
    end interface

    h = impl(comm)
    if (present(ierr)) then
      if (c_associated(h)) then
        ierr = 0  ! 0 == ESIO_SUCCESS
      else
        ierr = 5  ! 5 == ESIO_EFAILED
      end if
    end if

  end subroutine esio_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_finalize (h, ierr)
    type(esio_handle), intent(in)            :: h
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    interface
      function impl (h) bind (C, name="esio_finalize")
        import
        integer(c_int)                       :: impl
        type(esio_handle), intent(in), value :: h
      end function impl
    end interface

    stat = impl(h)
    if (present(ierr)) ierr = stat

  end subroutine esio_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_layout_count (count, ierr)

    integer, intent(out)           :: count
    integer, intent(out), optional :: ierr

    interface
      function impl () bind (C, name="esio_layout_count")
        import
        integer(c_int) :: impl
      end function impl
    end interface

    count = impl()
    if (present(ierr)) ierr = 0

  end subroutine esio_layout_count

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_layout_get (h, layout_index, ierr)

    type(esio_handle), intent(in)            :: h
    integer,           intent(out)           :: layout_index
    integer,           intent(out), optional :: ierr

    interface
      function impl (h) bind (C, name="esio_layout_get")
        import
        integer(c_int)                       :: impl
        type(esio_handle), intent(in), value :: h
      end function impl
    end interface

    layout_index = impl(h)
    if (present(ierr)) ierr = 0  ! FIXME: See bug #1178

  end subroutine esio_layout_get

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_layout_set (h, layout_index, ierr)

    type(esio_handle), intent(in)            :: h
    integer,           intent(in)            :: layout_index
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    interface
      function impl (h, layout_index) bind (C, name="esio_layout_set")
        import
        integer(c_int)                       :: impl
        type(esio_handle), intent(in), value :: h
        integer(c_int),    intent(in), value :: layout_index
      end function impl
    end interface

    stat = impl(h, layout_index)
    if (present(ierr)) ierr = stat

  end subroutine esio_layout_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_file_create (h, file, overwrite, ierr)

    type(esio_handle), intent(in)            :: h
    character(len=*),  intent(in)            :: file
    logical,           intent(in)            :: overwrite
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    interface
      function impl (h, file, overwrite) bind (C, name="esio_file_create")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: h
        character(len=1,kind=c_char), intent(in)        :: file(*)
        integer(c_int),               intent(in), value :: overwrite
      end function impl
    end interface

    stat = impl(h, f_c_string(file), f_c_logical(overwrite))
    if (present(ierr)) ierr = stat

  end subroutine esio_file_create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_file_open (h, file, readwrite, ierr)

    type(esio_handle), intent(in)            :: h
    character(len=*),  intent(in)            :: file
    logical,           intent(in)            :: readwrite
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    interface
      function impl (h, file, readwrite) bind (C, name="esio_file_open")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: h
        character(len=1,kind=c_char), intent(in)        :: file(*)
        integer(c_int),               intent(in), value :: readwrite
      end function impl
    end interface

    stat = impl(h, f_c_string(file), f_c_logical(readwrite))
    if (present(ierr)) ierr = stat

  end subroutine esio_file_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_file_flush (h, ierr)

    type(esio_handle), intent(in)            :: h
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    interface
      function impl (h) bind (C, name="esio_file_flush")
        import
        integer(c_int)                       :: impl
        type(esio_handle), intent(in), value :: h
      end function impl
    end interface

    stat = impl(h)
    if (present(ierr)) ierr = stat

  end subroutine esio_file_flush

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_file_close (h, ierr)

    type(esio_handle), intent(in)            :: h
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    interface
      function impl (h) bind (C, name="esio_file_close")
        import
        integer(c_int)                       :: impl
        type(esio_handle), intent(in), value :: h
      end function impl
    end interface

    stat = impl(h)
    if (present(ierr)) ierr = stat

  end subroutine esio_file_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_field_size (h, name, aglobal, bglobal, cglobal, ierr)

    type(esio_handle), intent(in)            :: h
    character(len=*),  intent(in)            :: name
    integer,           intent(out)           :: aglobal
    integer,           intent(out)           :: bglobal
    integer,           intent(out)           :: cglobal
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    integer(c_int) :: tmp_c, tmp_b, tmp_a

    interface
      function impl (h, name, cglobal, bglobal, aglobal)  &
                    bind (C, name="esio_field_size")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: h
        character(len=1,kind=c_char), intent(in)        :: name(*)
        integer(c_int),               intent(inout)     :: cglobal
        integer(c_int),               intent(inout)     :: bglobal
        integer(c_int),               intent(inout)     :: aglobal
      end function impl
    end interface

!   Note reordering Fortran's (a, b, c) to C's (c, b, a)
    stat = impl(h, f_c_string(name), tmp_c, tmp_b, tmp_a)
    aglobal = tmp_a
    bglobal = tmp_b
    cglobal = tmp_c
    if (present(ierr)) ierr = stat

  end subroutine esio_field_size


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_field_write_double (h, name, field,                   &
                                      aglobal, astart, alocal, astride, &
                                      bglobal, bstart, blocal, bstride, &
                                      cglobal, cstart, clocal, cstride, &
                                      ierr)

    type(esio_handle), intent(in) :: h
    character(len=*),  intent(in) :: name
    real(c_double),    intent(in) :: field(*)
    integer,           intent(in) :: aglobal, astart, alocal, astride
    integer,           intent(in) :: bglobal, bstart, blocal, bstride
    integer,           intent(in) :: cglobal, cstart, clocal, cstride
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    interface
      function impl (h, name, field,                   &
                     cglobal, cstart, clocal, cstride, &
                     bglobal, bstart, blocal, bstride, &
                     aglobal, astart, alocal, astride) &
                     bind (C, name="esio_field_write_double")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: h
        character(len=1,kind=c_char), intent(in)        :: name(*)
        real(c_double),               intent(in)        :: field(*)
        integer(c_int),               intent(in), value :: cglobal, &
                                                           cstart,  &
                                                           clocal,  &
                                                           cstride
        integer(c_int),               intent(in), value :: bglobal, &
                                                           bstart,  &
                                                           blocal,  &
                                                           bstride
        integer(c_int),               intent(in), value :: aglobal, &
                                                           astart,  &
                                                           alocal,  &
                                                           astride
      end function impl
    end interface

!   Note conversion from one- to zero-based starting offsets
!   Note reordering Fortran's (a, b, c) to C's (c, b, a)
    stat = impl(h, f_c_string(name), field,           &
                cglobal, cstart - 1, clocal, cstride, &
                bglobal, bstart - 1, blocal, bstride, &
                aglobal, astart - 1, alocal, astride)
    if (present(ierr)) ierr = stat

  end subroutine esio_field_write_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_field_write_single (h, name, field,                   &
                                      aglobal, astart, alocal, astride, &
                                      bglobal, bstart, blocal, bstride, &
                                      cglobal, cstart, clocal, cstride, &
                                      ierr)

    type(esio_handle), intent(in) :: h
    character(len=*),  intent(in) :: name
    real(c_float),     intent(in) :: field(*)
    integer,           intent(in) :: aglobal, astart, alocal, astride
    integer,           intent(in) :: bglobal, bstart, blocal, bstride
    integer,           intent(in) :: cglobal, cstart, clocal, cstride
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    interface
      function impl (h, name, field,                   &
                     cglobal, cstart, clocal, cstride, &
                     bglobal, bstart, blocal, bstride, &
                     aglobal, astart, alocal, astride) &
                     bind (C, name="esio_field_write_float")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: h
        character(len=1,kind=c_char), intent(in)        :: name(*)
        real(c_float),                intent(in)        :: field(*)
        integer(c_int),               intent(in), value :: cglobal, &
                                                           cstart,  &
                                                           clocal,  &
                                                           cstride
        integer(c_int),               intent(in), value :: bglobal, &
                                                           bstart,  &
                                                           blocal,  &
                                                           bstride
        integer(c_int),               intent(in), value :: aglobal, &
                                                           astart,  &
                                                           alocal,  &
                                                           astride
      end function impl
    end interface

!   Note conversion from one- to zero-based starting offsets
!   Note reordering Fortran's (a, b, c) to C's (c, b, a)
    stat = impl(h, f_c_string(name), field,           &
                cglobal, cstart - 1, clocal, cstride, &
                bglobal, bstart - 1, blocal, bstride, &
                aglobal, astart - 1, alocal, astride)
    if (present(ierr)) ierr = stat

  end subroutine esio_field_write_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_field_read_double (h, name, field,                   &
                                     aglobal, astart, alocal, astride, &
                                     bglobal, bstart, blocal, bstride, &
                                     cglobal, cstart, clocal, cstride, &
                                     ierr)

    type(esio_handle), intent(in)  :: h
    character(len=*),  intent(in)  :: name
    real(c_double),    intent(out) :: field(*)
    integer,           intent(in)  :: aglobal, astart, alocal, astride
    integer,           intent(in)  :: bglobal, bstart, blocal, bstride
    integer,           intent(in)  :: cglobal, cstart, clocal, cstride
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    interface
      function impl (h, name, field,                   &
                     cglobal, cstart, clocal, cstride, &
                     bglobal, bstart, blocal, bstride, &
                     aglobal, astart, alocal, astride) &
                     bind (C, name="esio_field_read_double")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: h
        character(len=1,kind=c_char), intent(in)        :: name(*)
        real(c_double),               intent(out)       :: field(*)
        integer(c_int),               intent(in), value :: cglobal, &
                                                           cstart,  &
                                                           clocal,  &
                                                           cstride
        integer(c_int),               intent(in), value :: bglobal, &
                                                           bstart,  &
                                                           blocal,  &
                                                           bstride
        integer(c_int),               intent(in), value :: aglobal, &
                                                           astart,  &
                                                           alocal,  &
                                                           astride
      end function impl
    end interface

!   Note conversion from one- to zero-based starting offsets
!   Note reordering Fortran's (a, b, c) to C's (c, b, a)
    stat = impl(h, f_c_string(name), field,           &
                cglobal, cstart - 1, clocal, cstride, &
                bglobal, bstart - 1, blocal, bstride, &
                aglobal, astart - 1, alocal, astride)
    if (present(ierr)) ierr = stat

  end subroutine esio_field_read_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_field_read_single (h, name, field,                    &
                                     aglobal, astart, alocal, astride,  &
                                     bglobal, bstart, blocal, bstride,  &
                                     cglobal, cstart, clocal, cstride,  &
                                     ierr)

    type(esio_handle), intent(in)  :: h
    character(len=*),  intent(in)  :: name
    real(c_float),     intent(out) :: field(*)
    integer,           intent(in)  :: aglobal, astart, alocal, astride
    integer,           intent(in)  :: bglobal, bstart, blocal, bstride
    integer,           intent(in)  :: cglobal, cstart, clocal, cstride
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    interface
      function impl (h, name, field,                   &
                     cglobal, cstart, clocal, cstride, &
                     bglobal, bstart, blocal, bstride, &
                     aglobal, astart, alocal, astride) &
                     bind (C, name="esio_field_read_float")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: h
        character(len=1,kind=c_char), intent(in)        :: name(*)
        real(c_float),                intent(out)       :: field(*)
        integer(c_int),               intent(in), value :: cglobal, &
                                                           cstart,  &
                                                           clocal,  &
                                                           cstride
        integer(c_int),               intent(in), value :: bglobal, &
                                                           bstart,  &
                                                           blocal,  &
                                                           bstride
        integer(c_int),               intent(in), value :: aglobal, &
                                                           astart,  &
                                                           alocal,  &
                                                           astride
      end function impl
    end interface

!   Note conversion from one- to zero-based starting offsets
!   Note reordering Fortran's (a, b, c) to C's (c, b, a)
    stat = impl(h, f_c_string(name), field,           &
                cglobal, cstart - 1, clocal, cstride, &
                bglobal, bstart - 1, blocal, bstride, &
                aglobal, astart - 1, alocal, astride)
    if (present(ierr)) ierr = stat

  end subroutine esio_field_read_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>Convert Fortran-style string to a C-style one.
!!Trailing whitespace is trimmed.
!!See http://fortranwiki.org/fortran/show/Generating+C+Interfaces#strings
  pure function f_c_string (f_string) result (c_string)

    use, intrinsic :: iso_c_binding, only: c_char, c_null_char

    implicit none

    character(len=*), intent(in) :: f_string
    character(len=1,kind=c_char) :: c_string(len_trim(f_string)+1)
    integer                      :: n, i

    n = len_trim(f_string)
    do i = 1, n
!     Trivial slice required; f_string(i) parsed as an invalid function call
      c_string(i) = f_string(i:i)
    end do
    c_string(n + 1) = c_null_char

  end function f_c_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>Convert Fortran \c logical to a C-style true/false \c int.
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
