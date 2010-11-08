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
                                         esio_handle => c_ptr

  implicit none  ! Nothing implicit

! C interoperation details are kept hidden from the client...
  private :: c_char, c_double, c_double_complex
  private :: c_float, c_float_complex, c_int, c_null_char
  private :: c_associated
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

! TODO Allow Fortran to use customizable error handling
! Error handling routine
  private :: esio_error

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

  subroutine esio_string_set (h, name, value, ierr)

    type(esio_handle), intent(in)            :: h
    character(len=*),  intent(in)            :: name
    character(len=*),  intent(in)            :: value
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    interface
      function impl (h, name, value) bind (C, name="esio_string_set")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: h
        character(len=1,kind=c_char), intent(in)        :: name(*)
        character(len=1,kind=c_char), intent(in)        :: value(*)
      end function impl
    end interface

    stat = impl(h, f_c_string(name), f_c_string(value))
    if (present(ierr)) ierr = stat

  end subroutine esio_string_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_string_get (h, name, value, ierr)

    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_f_pointer

    type(esio_handle), intent(in)            :: h
    character(len=*),  intent(in)            :: name
    character(len=*),  intent(out)           :: value
    integer,           intent(out), optional :: ierr

    type(c_ptr)                              :: tmp_p
    character(len=1,kind=c_char), pointer    :: tmp_str(:)
    integer                                  :: i, n(1)

!   The C implementation which returns newly allocated memory
    interface
      function impl (h, name) bind (C, name="esio_string_get")
        import
        type(c_ptr)                                     :: impl
        type(esio_handle),            intent(in), value :: h
        character(len=1,kind=c_char), intent(in)        :: name(*)
      end function impl
    end interface

!   An (unavoidable?) abuse of the ISO_C_BINDING to get at C's free
    interface
      subroutine c_free (p) bind (C, name="free")
        import
        type(c_ptr), intent(in) :: p
      end subroutine c_free
    end interface

    tmp_p = impl(h, f_c_string(name))
    if (c_associated(tmp_p)) then
      value = ''
      n(1) = len(value)
      call c_f_pointer(tmp_p, tmp_str, n)
      do i = 1, n(1)
        if (tmp_str(i) == c_null_char) exit
        value(i:i) = tmp_str(i)
      end do
      call c_free(tmp_p)
      if (present(ierr)) ierr = 0  ! 0 == ESIO_SUCCESS
    else
      if (present(ierr)) then
        ierr = 5  ! 5 == ESIO_EFAILED
      else
        call esio_error('esio_string_get failed but ierr was not supplied', &
                        __FILE__, __LINE__, 5)
        call abort
      endif
    endif

  end subroutine esio_string_get

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_attribute_write_double
#define FINTENT intent(in)
#define CTYPE real(c_double)
#define CBINDNAME "esio_attribute_write_double"
#include "attribute.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_attribute_write_single
#define FINTENT intent(in)
#define CTYPE real(c_float)
#define CBINDNAME "esio_attribute_write_float"
#include "attribute.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_attribute_write_integer
#define FINTENT intent(in)
#define CTYPE integer(c_int)
#define CBINDNAME "esio_attribute_write_int"
#include "attribute.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_attribute_read_double
#define FINTENT intent(out)
#define CTYPE real(c_double)
#define CBINDNAME "esio_attribute_read_double"
#include "attribute.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_attribute_read_single
#define FINTENT intent(out)
#define CTYPE real(c_float)
#define CBINDNAME "esio_attribute_read_float"
#include "attribute.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_attribute_read_integer
#define FINTENT intent(out)
#define CTYPE integer(c_int)
#define CBINDNAME "esio_attribute_read_int"
#include "attribute.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! No esio_attribute_size in C interface so none in Fortran either

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_attribute_writev_double
#define FINTENT intent(in)
#define CTYPE real(c_double)
#define CBINDNAME "esio_attribute_writev_double"
#include "attribute.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_attribute_writev_single
#define FINTENT intent(in)
#define CTYPE real(c_float)
#define CBINDNAME "esio_attribute_writev_float"
#include "attribute.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_attribute_writev_integer
#define FINTENT intent(in)
#define CTYPE integer(c_int)
#define CBINDNAME "esio_attribute_writev_int"
#include "attribute.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_attribute_readv_double
#define FINTENT intent(out)
#define CTYPE real(c_double)
#define CBINDNAME "esio_attribute_readv_double"
#include "attribute.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_attribute_readv_single
#define FINTENT intent(out)
#define CTYPE real(c_float)
#define CBINDNAME "esio_attribute_readv_float"
#include "attribute.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_attribute_readv_integer
#define FINTENT intent(out)
#define CTYPE integer(c_int)
#define CBINDNAME "esio_attribute_readv_int"
#include "attribute.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_attribute_sizev (h, name, ncomponents, ierr)

    type(esio_handle), intent(in)            :: h
    character(len=*),  intent(in)            :: name
    integer,           intent(out)           :: ncomponents
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    integer(c_int) :: tmp_ncomponents

    interface
      function impl (h, name, ncomponents)  &
                    bind (C, name="esio_attribute_sizev")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: h
        character(len=1,kind=c_char), intent(in)        :: name(*)
        integer(c_int),               intent(inout)     :: ncomponents
      end function impl
    end interface

    stat = impl(h, f_c_string(name), tmp_ncomponents)
    ncomponents = tmp_ncomponents
    if (present(ierr)) ierr = stat

  end subroutine esio_attribute_sizev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_line_write_double
#define FINTENT intent(in)
#define CTYPE real(c_double)
#define CBINDNAME "esio_line_write_double"
#include "line.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_line_write_single
#define FINTENT intent(in)
#define CTYPE real(c_float)
#define CBINDNAME "esio_line_write_float"
#include "line.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_line_write_integer
#define FINTENT intent(in)
#define CTYPE integer(c_int)
#define CBINDNAME "esio_line_write_int"
#include "line.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_line_read_double
#define FINTENT intent(out)
#define CTYPE real(c_double)
#define CBINDNAME "esio_line_read_double"
#include "line.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_line_read_single
#define FINTENT intent(out)
#define CTYPE real(c_float)
#define CBINDNAME "esio_line_read_float"
#include "line.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_line_read_integer
#define FINTENT intent(out)
#define CTYPE integer(c_int)
#define CBINDNAME "esio_line_read_int"
#include "line.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_line_size (h, name, aglobal, ierr)

    type(esio_handle), intent(in)            :: h
    character(len=*),  intent(in)            :: name
    integer,           intent(out)           :: aglobal
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    integer(c_int) :: tmp_a

    interface
      function impl (h, name, aglobal) bind (C, name="esio_line_size")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: h
        character(len=1,kind=c_char), intent(in)        :: name(*)
        integer(c_int),               intent(inout)     :: aglobal
      end function impl
    end interface

    stat = impl(h, f_c_string(name), tmp_a)
    aglobal = tmp_a
    if (present(ierr)) ierr = stat

  end subroutine esio_line_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_line_writev_double
#define FINTENT intent(in)
#define CTYPE real(c_double)
#define CBINDNAME "esio_line_writev_double"
#include "line.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_line_writev_single
#define FINTENT intent(in)
#define CTYPE real(c_float)
#define CBINDNAME "esio_line_writev_float"
#include "line.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_line_writev_integer
#define FINTENT intent(in)
#define CTYPE integer(c_int)
#define CBINDNAME "esio_line_writev_int"
#include "line.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_line_readv_double
#define FINTENT intent(out)
#define CTYPE real(c_double)
#define CBINDNAME "esio_line_readv_double"
#include "line.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_line_readv_single
#define FINTENT intent(out)
#define CTYPE real(c_float)
#define CBINDNAME "esio_line_readv_float"
#include "line.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_line_readv_integer
#define FINTENT intent(out)
#define CTYPE integer(c_int)
#define CBINDNAME "esio_line_readv_int"
#include "line.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_line_sizev (h, name, aglobal, ncomponents, ierr)

    type(esio_handle), intent(in)            :: h
    character(len=*),  intent(in)            :: name
    integer,           intent(out)           :: aglobal
    integer,           intent(out)           :: ncomponents
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    integer(c_int) :: tmp_a, tmp_ncomponents

    interface
      function impl (h, name, aglobal, ncomponents)  &
                    bind (C, name="esio_line_sizev")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: h
        character(len=1,kind=c_char), intent(in)        :: name(*)
        integer(c_int),               intent(inout)     :: aglobal
        integer(c_int),               intent(inout)     :: ncomponents
      end function impl
    end interface

    stat = impl(h, f_c_string(name), tmp_a, tmp_ncomponents)
    aglobal = tmp_a
    ncomponents = tmp_ncomponents
    if (present(ierr)) ierr = stat

  end subroutine esio_line_sizev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_plane_write_double
#define FINTENT intent(in)
#define CTYPE real(c_double)
#define CBINDNAME "esio_plane_write_double"
#include "plane.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_plane_write_single
#define FINTENT intent(in)
#define CTYPE real(c_float)
#define CBINDNAME "esio_plane_write_float"
#include "plane.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_plane_write_integer
#define FINTENT intent(in)
#define CTYPE integer(c_int)
#define CBINDNAME "esio_plane_write_int"
#include "plane.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_plane_read_double
#define FINTENT intent(out)
#define CTYPE real(c_double)
#define CBINDNAME "esio_plane_read_double"
#include "plane.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_plane_read_single
#define FINTENT intent(out)
#define CTYPE real(c_float)
#define CBINDNAME "esio_plane_read_float"
#include "plane.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_plane_read_integer
#define FINTENT intent(out)
#define CTYPE integer(c_int)
#define CBINDNAME "esio_plane_read_int"
#include "plane.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_plane_size (h, name, aglobal, bglobal, ierr)

    type(esio_handle), intent(in)            :: h
    character(len=*),  intent(in)            :: name
    integer,           intent(out)           :: aglobal
    integer,           intent(out)           :: bglobal
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    integer(c_int) :: tmp_b, tmp_a

    interface
      function impl (h, name, bglobal, aglobal)  &
                    bind (C, name="esio_plane_size")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: h
        character(len=1,kind=c_char), intent(in)        :: name(*)
        integer(c_int),               intent(inout)     :: bglobal
        integer(c_int),               intent(inout)     :: aglobal
      end function impl
    end interface

!   Note reordering Fortran's (a, b) to C's (b, a)
    stat = impl(h, f_c_string(name), tmp_b, tmp_a)
    aglobal = tmp_a
    bglobal = tmp_b
    if (present(ierr)) ierr = stat

  end subroutine esio_plane_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_plane_writev_double
#define FINTENT intent(in)
#define CTYPE real(c_double)
#define CBINDNAME "esio_plane_writev_double"
#include "plane.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_plane_writev_single
#define FINTENT intent(in)
#define CTYPE real(c_float)
#define CBINDNAME "esio_plane_writev_float"
#include "plane.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_plane_writev_integer
#define FINTENT intent(in)
#define CTYPE integer(c_int)
#define CBINDNAME "esio_plane_writev_int"
#include "plane.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_plane_readv_double
#define FINTENT intent(out)
#define CTYPE real(c_double)
#define CBINDNAME "esio_plane_readv_double"
#include "plane.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_plane_readv_single
#define FINTENT intent(out)
#define CTYPE real(c_float)
#define CBINDNAME "esio_plane_readv_float"
#include "plane.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_plane_readv_integer
#define FINTENT intent(out)
#define CTYPE integer(c_int)
#define CBINDNAME "esio_plane_readv_int"
#include "plane.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_plane_sizev (h, name, aglobal, bglobal, ncomponents, ierr)

    type(esio_handle), intent(in)            :: h
    character(len=*),  intent(in)            :: name
    integer,           intent(out)           :: aglobal
    integer,           intent(out)           :: bglobal
    integer,           intent(out)           :: ncomponents
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    integer(c_int) :: tmp_b, tmp_a, tmp_ncomponents

    interface
      function impl (h, name, bglobal, aglobal, ncomponents)  &
                    bind (C, name="esio_plane_sizev")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: h
        character(len=1,kind=c_char), intent(in)        :: name(*)
        integer(c_int),               intent(inout)     :: bglobal
        integer(c_int),               intent(inout)     :: aglobal
        integer(c_int),               intent(inout)     :: ncomponents
      end function impl
    end interface

!   Note reordering Fortran's (a, b) to C's (b, a)
    stat = impl(h, f_c_string(name), tmp_b, tmp_a, tmp_ncomponents)
    aglobal = tmp_a
    bglobal = tmp_b
    ncomponents = tmp_ncomponents
    if (present(ierr)) ierr = stat

  end subroutine esio_plane_sizev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_field_write_double
#define FINTENT intent(in)
#define CTYPE real(c_double)
#define CBINDNAME "esio_field_write_double"
#include "field.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_field_write_single
#define FINTENT intent(in)
#define CTYPE real(c_float)
#define CBINDNAME "esio_field_write_float"
#include "field.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_field_write_integer
#define FINTENT intent(in)
#define CTYPE integer(c_int)
#define CBINDNAME "esio_field_write_int"
#include "field.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_field_read_double
#define FINTENT intent(out)
#define CTYPE real(c_double)
#define CBINDNAME "esio_field_read_double"
#include "field.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_field_read_single
#define FINTENT intent(out)
#define CTYPE real(c_float)
#define CBINDNAME "esio_field_read_float"
#include "field.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define FNAME esio_field_read_integer
#define FINTENT intent(out)
#define CTYPE integer(c_int)
#define CBINDNAME "esio_field_read_int"
#include "field.f90"

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

#define VECTORVALUED
#define FNAME esio_field_writev_double
#define FINTENT intent(in)
#define CTYPE real(c_double)
#define CBINDNAME "esio_field_writev_double"
#include "field.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_field_writev_single
#define FINTENT intent(in)
#define CTYPE real(c_float)
#define CBINDNAME "esio_field_writev_float"
#include "field.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_field_writev_integer
#define FINTENT intent(in)
#define CTYPE integer(c_int)
#define CBINDNAME "esio_field_writev_int"
#include "field.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_field_readv_double
#define FINTENT intent(out)
#define CTYPE real(c_double)
#define CBINDNAME "esio_field_readv_double"
#include "field.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_field_readv_single
#define FINTENT intent(out)
#define CTYPE real(c_float)
#define CBINDNAME "esio_field_readv_float"
#include "field.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VECTORVALUED
#define FNAME esio_field_readv_integer
#define FINTENT intent(out)
#define CTYPE integer(c_int)
#define CBINDNAME "esio_field_readv_int"
#include "field.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_field_sizev (h, name, aglobal, bglobal, cglobal, &
                               ncomponents, ierr)

    type(esio_handle), intent(in)            :: h
    character(len=*),  intent(in)            :: name
    integer,           intent(out)           :: aglobal
    integer,           intent(out)           :: bglobal
    integer,           intent(out)           :: cglobal
    integer,           intent(out)           :: ncomponents
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    integer(c_int) :: tmp_c, tmp_b, tmp_a, tmp_ncomponents

    interface
      function impl (h, name, cglobal, bglobal, aglobal, ncomponents)  &
                    bind (C, name="esio_field_sizev")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: h
        character(len=1,kind=c_char), intent(in)        :: name(*)
        integer(c_int),               intent(inout)     :: cglobal
        integer(c_int),               intent(inout)     :: bglobal
        integer(c_int),               intent(inout)     :: aglobal
        integer(c_int),               intent(inout)     :: ncomponents
      end function impl
    end interface

!   Note reordering Fortran's (a, b, c) to C's (c, b, a)
    stat = impl(h, f_c_string(name), tmp_c, tmp_b, tmp_a, tmp_ncomponents)
    aglobal = tmp_a
    bglobal = tmp_b
    cglobal = tmp_c
    ncomponents = tmp_ncomponents
    if (present(ierr)) ierr = stat

  end subroutine esio_field_sizev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_error (reason, file, line, esio_errno)

    character(len=*),  intent(in) :: reason
    character(len=*),  intent(in) :: file
    integer,           intent(in) :: line
    integer,           intent(in) :: esio_errno

    interface
      subroutine c_impl (reason, file, line, esio_errno)  &
                 bind (C, name="esio_error")
        import
        character(len=1,kind=c_char), intent(in)        :: reason(*)
        character(len=1,kind=c_char), intent(in)        :: file(*)
        integer(c_int),               intent(in), value :: line
        integer(c_int),               intent(in), value :: esio_errno
      end subroutine
    end interface

    call c_impl(f_c_string(reason), f_c_string(file), line, esio_errno)

  end subroutine esio_error

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
