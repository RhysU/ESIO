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

!> \file
!! Provides ESIO's Fortran-based public API following the library's
!! \ref concepts "usage concepts".  All methods in this header
!! invoke ESIO's \ref error.h "error handling mechanisms" on failure.

module esio

! Use ISO_C_BINDING for interoperation with ESIO's C implementation.
! Rename C_PTR to ESIO_STATE to make the handle object types more opaque.
! State traverses the language boundary as TYPE(ESIO_STATE)s.
! Renaming also prevents collision if caller uses ISO_C_BINDING directly.
  use, intrinsic :: iso_c_binding, only: c_char,             &
                                         c_double,           &
                                         c_double_complex,   &
                                         c_float,            &
                                         c_float_complex,    &
                                         c_int,              &
                                         c_null_char,        &
                                         esio_handle => c_ptr

  implicit none  ! Nothing implicit
  private        ! Everything default private

! C interoperation details are kept hidden from the client...
! ... with the exception of our opaque handle object type.
  public :: esio_handle

! Public API
  public :: esio_initialize, esio_finalize
  public :: esio_layout_count, esio_layout_get, esio_layout_set
  public :: esio_file_create, esio_file_open, esio_file_close
  public :: esio_field_size
  public :: esio_field_write_double, esio_field_write_single
  public :: esio_field_read_double, esio_field_read_single

! Generic, precision-agnostic interfaces atop the public API
  public :: esio_field_write
  interface esio_field_write
    module procedure esio_field_write_double
    module procedure esio_field_write_single
  end interface

  public :: esio_field_read
  interface esio_field_read
    module procedure esio_field_read_double
    module procedure esio_field_read_single
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(esio_handle) function esio_initialize (comm)

    integer, intent(in) :: comm

    ! See C routine esio_initialize_fortran re: MPI communicator interoperation
    interface
      function impl (comm) bind (C, name="esio_initialize_fortran")
        import
        type(esio_handle)          :: impl
        integer, intent(in), value :: comm  ! Note integer not integer(c_int)
      end function impl
    end interface

    esio_initialize = impl(comm)

  end function esio_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_finalize (handle)

    type(esio_handle), intent(in) :: handle

    interface
      function impl (handle) bind (C, name="esio_finalize")
        import
        integer(c_int)                       :: impl
        type(esio_handle), intent(in), value :: handle
      end function impl
    end interface

    esio_finalize = impl(handle)

  end function esio_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_layout_count ()

    interface
      function impl () bind (C, name="esio_layout_count")
        import
        integer(c_int) :: impl
      end function impl
    end interface

    esio_layout_count = impl()

  end function esio_layout_count

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_layout_get (handle)

    type(esio_handle), intent(in) :: handle

    interface
      function impl (handle) bind (C, name="esio_layout_get")
        import
        integer(c_int)                       :: impl
        type(esio_handle), intent(in), value :: handle
      end function impl
    end interface

    esio_layout_get = impl(handle)

  end function esio_layout_get

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_layout_set (handle, layout_tag)

    type(esio_handle), intent(in) :: handle
    integer,           intent(in) :: layout_tag

    interface
      function impl (handle, layout_tag) bind (C, name="esio_layout_set")
        import
        integer(c_int)                       :: impl
        type(esio_handle), intent(in), value :: handle
        integer(c_int),    intent(in), value :: layout_tag
      end function impl
    end interface

    esio_layout_set = impl(handle, layout_tag)

  end function esio_layout_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_file_create (handle, file, overwrite)

    type(esio_handle), intent(in) :: handle
    character(len=*),  intent(in) :: file
    logical,           intent(in) :: overwrite

    interface
      function impl (handle, file, overwrite) bind (C, name="esio_file_create")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: handle
        character(len=1,kind=c_char), intent(in)        :: file(*)
        integer(c_int),               intent(in), value :: overwrite
      end function impl
    end interface

    esio_file_create = impl(handle, f_c_string(file), f_c_logical(overwrite))

  end function esio_file_create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_file_open (handle, file, readwrite)

    type(esio_handle), intent(in) :: handle
    character(len=*),  intent(in) :: file
    logical,           intent(in) :: readwrite

    interface
      function impl (handle, file, readwrite) bind (C, name="esio_file_open")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: handle
        character(len=1,kind=c_char), intent(in)        :: file(*)
        integer(c_int),               intent(in), value :: readwrite
      end function impl
    end interface

    esio_file_open = impl(handle, f_c_string(file), f_c_logical(readwrite))

  end function esio_file_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_file_close (handle)

    type(esio_handle), intent(in) :: handle

    interface
      function impl (handle) bind (C, name="esio_file_close")
        import
        integer(c_int)                       :: impl
        type(esio_handle), intent(in), value :: handle
      end function impl
    end interface

    esio_file_close = impl(handle)

  end function esio_file_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_field_size (handle, name, aglobal, bglobal, cglobal)

    type(esio_handle), intent(in)  :: handle
    character(len=*),  intent(in)  :: name
    integer,           intent(out) :: aglobal
    integer,           intent(out) :: bglobal
    integer,           intent(out) :: cglobal

    integer(c_int) :: tmp_c, tmp_b, tmp_a

    interface
      function impl (handle, name, cglobal, bglobal, aglobal)  &
                    bind (C, name="esio_field_size")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: handle
        character(len=1,kind=c_char), intent(in)        :: name(*)
        integer(c_int),               intent(inout)     :: cglobal
        integer(c_int),               intent(inout)     :: bglobal
        integer(c_int),               intent(inout)     :: aglobal
      end function impl
    end interface

!   Note reordering Fortran's (a, b, c) to C's (c, b, a)
    esio_field_size = impl(handle, f_c_string(name), tmp_c, tmp_b, tmp_a)
    aglobal = tmp_a
    bglobal = tmp_b
    cglobal = tmp_c

  end function esio_field_size


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_field_write_double (handle, name, field,              &
                                            aglobal, astart, alocal, astride, &
                                            bglobal, bstart, blocal, bstride, &
                                            cglobal, cstart, clocal, cstride)

    type(esio_handle), intent(in) :: handle
    character(len=*),  intent(in) :: name
    real(c_double),    intent(in) :: field(*)
    integer,           intent(in) :: aglobal, astart, alocal, astride
    integer,           intent(in) :: bglobal, bstart, blocal, bstride
    integer,           intent(in) :: cglobal, cstart, clocal, cstride

    interface
      function impl (handle, name, field,              &
                     cglobal, cstart, clocal, cstride, &
                     bglobal, bstart, blocal, bstride, &
                     aglobal, astart, alocal, astride) &
                     bind (C, name="esio_field_write_double")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: handle
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
    esio_field_write_double = impl(handle, f_c_string(name), field,      &
                                   cglobal, cstart - 1, clocal, cstride, &
                                   bglobal, bstart - 1, blocal, bstride, &
                                   aglobal, astart - 1, alocal, astride)

  end function esio_field_write_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_field_write_single (handle, name, field,              &
                                            aglobal, astart, alocal, astride, &
                                            bglobal, bstart, blocal, bstride, &
                                            cglobal, cstart, clocal, cstride)

    type(esio_handle), intent(in) :: handle
    character(len=*),  intent(in) :: name
    real(c_float),     intent(in) :: field(*)
    integer,           intent(in) :: aglobal, astart, alocal, astride
    integer,           intent(in) :: bglobal, bstart, blocal, bstride
    integer,           intent(in) :: cglobal, cstart, clocal, cstride

    interface
      function impl (handle, name, field,              &
                     cglobal, cstart, clocal, cstride, &
                     bglobal, bstart, blocal, bstride, &
                     aglobal, astart, alocal, astride) &
                     bind (C, name="esio_field_write_float")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: handle
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
    esio_field_write_single = impl(handle, f_c_string(name), field,      &
                                   cglobal, cstart - 1, clocal, cstride, &
                                   bglobal, bstart - 1, blocal, bstride, &
                                   aglobal, astart - 1, alocal, astride)

  end function esio_field_write_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_field_read_double (handle, name, field,              &
                                           aglobal, astart, alocal, astride, &
                                           bglobal, bstart, blocal, bstride, &
                                           cglobal, cstart, clocal, cstride)

    type(esio_handle), intent(in)  :: handle
    character(len=*),  intent(in)  :: name
    real(c_double),    intent(out) :: field(*)
    integer,           intent(in)  :: aglobal, astart, alocal, astride
    integer,           intent(in)  :: bglobal, bstart, blocal, bstride
    integer,           intent(in)  :: cglobal, cstart, clocal, cstride

    interface
      function impl (handle, name, field,              &
                     cglobal, cstart, clocal, cstride, &
                     bglobal, bstart, blocal, bstride, &
                     aglobal, astart, alocal, astride) &
                     bind (C, name="esio_field_read_double")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: handle
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
    esio_field_read_double = impl(handle, f_c_string(name), field,      &
                                  cglobal, cstart - 1, clocal, cstride, &
                                  bglobal, bstart - 1, blocal, bstride, &
                                  aglobal, astart - 1, alocal, astride)

  end function esio_field_read_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_field_read_single (handle, name, field,               &
                                           aglobal, astart, alocal, astride,  &
                                           bglobal, bstart, blocal, bstride,  &
                                           cglobal, cstart, clocal, cstride)

    type(esio_handle), intent(in)  :: handle
    character(len=*),  intent(in)  :: name
    real(c_float),     intent(out) :: field(*)
    integer,           intent(in)  :: aglobal, astart, alocal, astride
    integer,           intent(in)  :: bglobal, bstart, blocal, bstride
    integer,           intent(in)  :: cglobal, cstart, clocal, cstride

    interface
      function impl (handle, name, field,              &
                     cglobal, cstart, clocal, cstride, &
                     bglobal, bstart, blocal, bstride, &
                     aglobal, astart, alocal, astride) &
                     bind (C, name="esio_field_read_float")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: handle
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
    esio_field_read_single = impl(handle, f_c_string(name), field,      &
                                  cglobal, cstart - 1, clocal, cstride, &
                                  bglobal, bstart - 1, blocal, bstride, &
                                  aglobal, astart - 1, alocal, astride)

  end function esio_field_read_single

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
