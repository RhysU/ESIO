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

! Generic, precision-agnostic attribute interfaces atop the public API
  interface esio_attribute_write
    module procedure esio_attribute_write_double
    module procedure esio_attribute_write_single
    module procedure esio_attribute_write_integer
  end interface

  interface esio_attribute_read
    module procedure esio_attribute_read_double
    module procedure esio_attribute_read_single
    module procedure esio_attribute_read_integer
  end interface

  interface esio_attribute_writev
    module procedure esio_attribute_writev_double
    module procedure esio_attribute_writev_single
    module procedure esio_attribute_writev_integer
  end interface

  interface esio_attribute_readv
    module procedure esio_attribute_readv_double
    module procedure esio_attribute_readv_single
    module procedure esio_attribute_readv_integer
  end interface

! Generic, precision-agnostic line interfaces atop the public API
  interface esio_line_write
    module procedure esio_line_write_double
    module procedure esio_line_write_single
    module procedure esio_line_write_integer
  end interface

  interface esio_line_read
    module procedure esio_line_read_double
    module procedure esio_line_read_single
    module procedure esio_line_read_integer
  end interface

  interface esio_line_writev
    module procedure esio_line_writev_double
    module procedure esio_line_writev_single
    module procedure esio_line_writev_integer
  end interface

  interface esio_line_readv
    module procedure esio_line_readv_double
    module procedure esio_line_readv_single
    module procedure esio_line_readv_integer
  end interface

! Generic, precision-agnostic plane interfaces atop the public API
  interface esio_plane_write
    module procedure esio_plane_write_double
    module procedure esio_plane_write_single
    module procedure esio_plane_write_integer
  end interface

  interface esio_plane_read
    module procedure esio_plane_read_double
    module procedure esio_plane_read_single
    module procedure esio_plane_read_integer
  end interface

  interface esio_plane_writev
    module procedure esio_plane_writev_double
    module procedure esio_plane_writev_single
    module procedure esio_plane_writev_integer
  end interface

  interface esio_plane_readv
    module procedure esio_plane_readv_double
    module procedure esio_plane_readv_single
    module procedure esio_plane_readv_integer
  end interface

! Generic, precision-agnostic field interfaces atop the public API
  interface esio_field_write
    module procedure esio_field_write_double
    module procedure esio_field_write_single
    module procedure esio_field_write_integer
  end interface

  interface esio_field_read
    module procedure esio_field_read_double
    module procedure esio_field_read_single
    module procedure esio_field_read_integer
  end interface

  interface esio_field_writev
    module procedure esio_field_writev_double
    module procedure esio_field_writev_single
    module procedure esio_field_writev_integer
  end interface

  interface esio_field_readv
    module procedure esio_field_readv_double
    module procedure esio_field_readv_single
    module procedure esio_field_readv_integer
  end interface


! TODO Allow Fortran to use customizable error handling
! Error handling routine
  private :: esio_error

! Internal helper routines
  private :: f_c_string, f_c_logical

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Initializing and finalizing a state handle
!! See \ref conceptshandles "handles" for the associated semantics.
!!@{

  subroutine esio_handle_initialize (h, comm, ierr)

    type(esio_handle), intent(out)           :: h
    integer,           intent(in)            :: comm
    integer,           intent(out), optional :: ierr

    ! See C routine esio_handle_initialize_fortran
    ! for details on MPI communicator interoperation
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (comm) bind (C, name="esio_handle_initialize_fortran")
        import
        type(esio_handle)          :: impl
        integer, intent(in), value :: comm  ! Note integer not integer(c_int)
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    h = impl(comm)
    if (present(ierr)) then
      if (c_associated(h)) then
        ierr = 0  ! 0 == ESIO_SUCCESS
      else
        ierr = 5  ! 5 == ESIO_EFAILED
      end if
    end if

  end subroutine esio_handle_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_handle_finalize (h, ierr)
    type(esio_handle), intent(in)            :: h
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (h) bind (C, name="esio_handle_finalize")
        import
        integer(c_int)                       :: impl
        type(esio_handle), intent(in), value :: h
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    stat = impl(h)
    if (present(ierr)) ierr = stat

  end subroutine esio_handle_finalize

!> @}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Opening and closing files
!! See \ref conceptsfiles "file concepts" for more details.
!!@{

  subroutine esio_file_create (h, file, overwrite, ierr)

    type(esio_handle), intent(in)            :: h
    character(len=*),  intent(in)            :: file
    logical,           intent(in)            :: overwrite
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (h, file, overwrite) bind (C, name="esio_file_create")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: h
        character(len=1,kind=c_char), intent(in)        :: file(*)
        integer(c_int),               intent(in), value :: overwrite
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

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

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (h, file, readwrite) bind (C, name="esio_file_open")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: h
        character(len=1,kind=c_char), intent(in)        :: file(*)
        integer(c_int),               intent(in), value :: readwrite
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    stat = impl(h, f_c_string(file), f_c_logical(readwrite))
    if (present(ierr)) ierr = stat

  end subroutine esio_file_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_file_flush (h, ierr)

    type(esio_handle), intent(in)            :: h
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (h) bind (C, name="esio_file_flush")
        import
        integer(c_int)                       :: impl
        type(esio_handle), intent(in), value :: h
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    stat = impl(h)
    if (present(ierr)) ierr = stat

  end subroutine esio_file_flush

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_file_close (h, ierr)

    type(esio_handle), intent(in)            :: h
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (h) bind (C, name="esio_file_close")
        import
        integer(c_int)                       :: impl
        type(esio_handle), intent(in), value :: h
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    stat = impl(h)
    if (present(ierr)) ierr = stat

  end subroutine esio_file_close

!!@}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Manipulating string-valued attributes
!! See \ref conceptsattributes "attributes concepts" for more details.
!!@{

  subroutine esio_string_set (h, name, value, ierr)

    type(esio_handle), intent(in)            :: h
    character(len=*),  intent(in)            :: name
    character(len=*),  intent(in)            :: value
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (h, name, value) bind (C, name="esio_string_set")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: h
        character(len=1,kind=c_char), intent(in)        :: name(*)
        character(len=1,kind=c_char), intent(in)        :: value(*)
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

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

!   The C implementation returns newly allocated memory
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (h, name) bind (C, name="esio_string_get")
        import
        type(c_ptr)                                     :: impl
        type(esio_handle),            intent(in), value :: h
        character(len=1,kind=c_char), intent(in)        :: name(*)
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

!   An (unavoidable?) abuse of the ISO_C_BINDING to get at C's free
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      subroutine c_free (p) bind (C, name="free")
        import
        type(c_ptr), intent(in) :: p
      end subroutine c_free
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

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

!!@}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Manipulating scalar-valued attributes
!! See \ref conceptsattributes "attribute concepts" for more details.
!!@{

subroutine esio_attribute_write_double (h, name, value, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(in)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_attribute_write_double"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_write_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_attribute_write_single (h, name, value, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(in)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_attribute_write_float"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_write_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_attribute_write_integer (h, name, value, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(in)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_attribute_write_int"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_write_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_attribute_read_double (h, name, value, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_attribute_read_double"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_read_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_attribute_read_single (h, name, value, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_attribute_read_float"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_read_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_attribute_read_integer (h, name, value, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_attribute_read_int"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_read_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! No esio_attribute_size in C interface so none in Fortran either

!!@}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Manipulating vector-valued attributes
!! See \ref conceptsattributes "attributes concepts" for more details.
!!@{

subroutine esio_attribute_writev_double (h, name, value, ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(in)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_attribute_writev_double"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_writev_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_attribute_writev_single (h, name, value, ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(in)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_attribute_writev_float"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_writev_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_attribute_writev_integer (h, name, value, ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(in)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_attribute_writev_int"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_writev_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_attribute_readv_double (h, name, value, ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_attribute_readv_double"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_readv_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_attribute_readv_single (h, name, value, ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_attribute_readv_float"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_readv_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_attribute_readv_integer (h, name, value, ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_attribute_readv_int"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_readv_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_attribute_sizev (h, name, ncomponents, ierr)

    type(esio_handle), intent(in)            :: h
    character(len=*),  intent(in)            :: name
    integer,           intent(out)           :: ncomponents
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    integer(c_int) :: tmp_ncomponents

#ifndef DOXYGEN_SHOULD_SKIP_THIS
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
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    stat = impl(h, f_c_string(name), tmp_ncomponents)
    ncomponents = tmp_ncomponents
    if (present(ierr)) ierr = stat

  end subroutine esio_attribute_sizev

!!@}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Manipulating distributed, one-dimensional, scalar-valued data
!! See \ref conceptslines "line concepts" for more details.
!!@{

subroutine esio_line_write_double (h, name, line,                    &
                                   aglobal, astart, alocal, astride, &
                                   ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(in)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_line_write_double"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_write_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_line_write_single (h, name, line,                    &
                                   aglobal, astart, alocal, astride, &
                                   ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(in)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_line_write_float"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_write_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_line_write_integer (h, name, line,                    &
                                    aglobal, astart, alocal, astride, &
                                    ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(in)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_line_write_int"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_write_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_line_read_double (h, name, line,                    &
                                  aglobal, astart, alocal, astride, &
                                  ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_line_read_double"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_read_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_line_read_single (h, name, line,                    &
                                  aglobal, astart, alocal, astride, &
                                  ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_line_read_float"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_read_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_line_read_integer (h, name, line,                    &
                                   aglobal, astart, alocal, astride, &
                                   ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_line_read_int"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_read_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_line_size (h, name, aglobal, ierr)

    type(esio_handle), intent(in)            :: h
    character(len=*),  intent(in)            :: name
    integer,           intent(out)           :: aglobal
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    integer(c_int) :: tmp_a

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (h, name, aglobal) bind (C, name="esio_line_size")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: h
        character(len=1,kind=c_char), intent(in)        :: name(*)
        integer(c_int),               intent(inout)     :: aglobal
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    stat = impl(h, f_c_string(name), tmp_a)
    aglobal = tmp_a
    if (present(ierr)) ierr = stat

  end subroutine esio_line_size

!!@}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Manipulating distributed, one-dimensional, vector-valued data
!! See \ref conceptslines "line concepts" for more details.
!!@{

subroutine esio_line_writev_double (h, name, line,                    &
                                    aglobal, astart, alocal, astride, &
                                    ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(in)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_line_writev_double"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_writev_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_line_writev_single (h, name, line,                    &
                                    aglobal, astart, alocal, astride, &
                                    ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(in)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_line_writev_float"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_writev_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_line_writev_integer (h, name, line,                    &
                                     aglobal, astart, alocal, astride, &
                                     ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(in)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_line_writev_int"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_writev_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_line_readv_double (h, name, line,                    &
                                   aglobal, astart, alocal, astride, &
                                   ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_line_readv_double"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_readv_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_line_readv_single (h, name, line,                    &
                                   aglobal, astart, alocal, astride, &
                                   ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_line_readv_float"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_readv_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_line_readv_integer (h, name, line,                    &
                                    aglobal, astart, alocal, astride, &
                                    ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_line_readv_int"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_readv_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_line_sizev (h, name, aglobal, ncomponents, ierr)

    type(esio_handle), intent(in)            :: h
    character(len=*),  intent(in)            :: name
    integer,           intent(out)           :: aglobal
    integer,           intent(out)           :: ncomponents
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    integer(c_int) :: tmp_a, tmp_ncomponents

#ifndef DOXYGEN_SHOULD_SKIP_THIS
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
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    stat = impl(h, f_c_string(name), tmp_a, tmp_ncomponents)
    aglobal = tmp_a
    ncomponents = tmp_ncomponents
    if (present(ierr)) ierr = stat

  end subroutine esio_line_sizev

!!@}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Manipulating distributed, two-dimensional, scalar-valued data
!! See \ref conceptsplanes "plane concepts" for more details.
!!@{

subroutine esio_plane_write_double (h, name, plane,                   &
                                    aglobal, astart, alocal, astride, &
                                    bglobal, bstart, blocal, bstride, &
                                    ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(in)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_plane_write_double"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_write_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_plane_write_single (h, name, plane,                   &
                                    aglobal, astart, alocal, astride, &
                                    bglobal, bstart, blocal, bstride, &
                                    ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(in)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_plane_write_float"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_write_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_plane_write_integer (h, name, plane,                   &
                                     aglobal, astart, alocal, astride, &
                                     bglobal, bstart, blocal, bstride, &
                                     ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(in)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_plane_write_int"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_write_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_plane_read_double (h, name, plane,                   &
                                   aglobal, astart, alocal, astride, &
                                   bglobal, bstart, blocal, bstride, &
                                   ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_plane_read_double"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_read_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_plane_read_single (h, name, plane,                   &
                                   aglobal, astart, alocal, astride, &
                                   bglobal, bstart, blocal, bstride, &
                                   ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_plane_read_float"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_read_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_plane_read_integer (h, name, plane,                   &
                                    aglobal, astart, alocal, astride, &
                                    bglobal, bstart, blocal, bstride, &
                                    ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_plane_read_int"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_read_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_plane_size (h, name, aglobal, bglobal, ierr)

    type(esio_handle), intent(in)            :: h
    character(len=*),  intent(in)            :: name
    integer,           intent(out)           :: aglobal
    integer,           intent(out)           :: bglobal
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    integer(c_int) :: tmp_b, tmp_a

#ifndef DOXYGEN_SHOULD_SKIP_THIS
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
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

!   Note reordering Fortran's (a, b) to C's (b, a)
    stat = impl(h, f_c_string(name), tmp_b, tmp_a)
    aglobal = tmp_a
    bglobal = tmp_b
    if (present(ierr)) ierr = stat

  end subroutine esio_plane_size

!!@}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Manipulating distributed, two-dimensional, vector-valued data
!! See \ref conceptsplanes "plane concepts" for more details.
!!@{

subroutine esio_plane_writev_double (h, name, plane,                   &
                                     aglobal, astart, alocal, astride, &
                                     bglobal, bstart, blocal, bstride, &
                                     ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(in)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_plane_writev_double"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_writev_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_plane_writev_single (h, name, plane,                   &
                                     aglobal, astart, alocal, astride, &
                                     bglobal, bstart, blocal, bstride, &
                                     ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(in)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_plane_writev_float"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_writev_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_plane_writev_integer (h, name, plane,                   &
                                      aglobal, astart, alocal, astride, &
                                      bglobal, bstart, blocal, bstride, &
                                      ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(in)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_plane_writev_int"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_writev_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_plane_readv_double (h, name, plane,                   &
                                    aglobal, astart, alocal, astride, &
                                    bglobal, bstart, blocal, bstride, &
                                    ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_plane_readv_double"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_readv_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_plane_readv_single (h, name, plane,                   &
                                    aglobal, astart, alocal, astride, &
                                    bglobal, bstart, blocal, bstride, &
                                    ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_plane_readv_float"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_readv_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_plane_readv_integer (h, name, plane,                   &
                                     aglobal, astart, alocal, astride, &
                                     bglobal, bstart, blocal, bstride, &
                                     ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_plane_readv_int"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_readv_integer

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

#ifndef DOXYGEN_SHOULD_SKIP_THIS
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
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

!   Note reordering Fortran's (a, b) to C's (b, a)
    stat = impl(h, f_c_string(name), tmp_b, tmp_a, tmp_ncomponents)
    aglobal = tmp_a
    bglobal = tmp_b
    ncomponents = tmp_ncomponents
    if (present(ierr)) ierr = stat

  end subroutine esio_plane_sizev

!!@}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Manipulating distributed, three-dimensional, scalar-valued data
!! See \ref conceptsfields "field concepts" for more details.
!!@{

subroutine esio_field_write_double (h, name, field,                   &
                                    aglobal, astart, alocal, astride, &
                                    bglobal, bstart, blocal, bstride, &
                                    cglobal, cstart, clocal, cstride, &
                                    ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(in)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_field_write_double"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_write_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_field_write_single (h, name, field,                   &
                                    aglobal, astart, alocal, astride, &
                                    bglobal, bstart, blocal, bstride, &
                                    cglobal, cstart, clocal, cstride, &
                                    ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(in)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_field_write_float"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_write_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_field_write_integer (h, name, field,                   &
                                     aglobal, astart, alocal, astride, &
                                     bglobal, bstart, blocal, bstride, &
                                     cglobal, cstart, clocal, cstride, &
                                     ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(in)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_field_write_int"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_write_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_field_read_double (h, name, field,                   &
                                   aglobal, astart, alocal, astride, &
                                   bglobal, bstart, blocal, bstride, &
                                   cglobal, cstart, clocal, cstride, &
                                   ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_field_read_double"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_read_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_field_read_single (h, name, field,                   &
                                   aglobal, astart, alocal, astride, &
                                   bglobal, bstart, blocal, bstride, &
                                   cglobal, cstart, clocal, cstride, &
                                   ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_field_read_float"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_read_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_field_read_integer (h, name, field,                   &
                                    aglobal, astart, alocal, astride, &
                                    bglobal, bstart, blocal, bstride, &
                                    cglobal, cstart, clocal, cstride, &
                                    ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_field_read_int"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_read_integer

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

#ifndef DOXYGEN_SHOULD_SKIP_THIS
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
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

!   Note reordering Fortran's (a, b, c) to C's (c, b, a)
    stat = impl(h, f_c_string(name), tmp_c, tmp_b, tmp_a)
    aglobal = tmp_a
    bglobal = tmp_b
    cglobal = tmp_c
    if (present(ierr)) ierr = stat

  end subroutine esio_field_size

!!@}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Manipulating distributed, three-dimensional, vector-valued data
!! See \ref conceptsfields "field concepts" for more details.
!!@{

subroutine esio_field_writev_double (h, name, field,                   &
                                     aglobal, astart, alocal, astride, &
                                     bglobal, bstart, blocal, bstride, &
                                     cglobal, cstart, clocal, cstride, &
                                     ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(in)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_field_writev_double"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_writev_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_field_writev_single (h, name, field,                   &
                                     aglobal, astart, alocal, astride, &
                                     bglobal, bstart, blocal, bstride, &
                                     cglobal, cstart, clocal, cstride, &
                                     ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(in)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_field_writev_float"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_writev_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_field_writev_integer (h, name, field,                   &
                                      aglobal, astart, alocal, astride, &
                                      bglobal, bstart, blocal, bstride, &
                                      cglobal, cstart, clocal, cstride, &
                                      ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(in)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_field_writev_int"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_writev_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_field_readv_double (h, name, field,                   &
                                    aglobal, astart, alocal, astride, &
                                    bglobal, bstart, blocal, bstride, &
                                    cglobal, cstart, clocal, cstride, &
                                    ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_field_readv_double"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_readv_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_field_readv_single (h, name, field,                   &
                                    aglobal, astart, alocal, astride, &
                                    bglobal, bstart, blocal, bstride, &
                                    cglobal, cstart, clocal, cstride, &
                                    ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_field_readv_float"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_readv_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_field_readv_integer (h, name, field,                   &
                                     aglobal, astart, alocal, astride, &
                                     bglobal, bstart, blocal, bstride, &
                                     cglobal, cstart, clocal, cstride, &
                                     ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_field_readv_int"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_readv_integer

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

#ifndef DOXYGEN_SHOULD_SKIP_THIS
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
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

!   Note reordering Fortran's (a, b, c) to C's (c, b, a)
    stat = impl(h, f_c_string(name), tmp_c, tmp_b, tmp_a, tmp_ncomponents)
    aglobal = tmp_a
    bglobal = tmp_b
    cglobal = tmp_c
    ncomponents = tmp_ncomponents
    if (present(ierr)) ierr = stat

  end subroutine esio_field_sizev

!!\@}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Querying and controlling field layout
!! See \ref conceptslayouts "layout concepts" for more details.
!!@{

  subroutine esio_field_layout_count (count, ierr)

    integer, intent(out)           :: count
    integer, intent(out), optional :: ierr

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl () bind (C, name="esio_field_layout_count")
        import
        integer(c_int) :: impl
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    count = impl()
    if (present(ierr)) ierr = 0

  end subroutine esio_field_layout_count

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_field_layout_get (h, layout_index, ierr)

    type(esio_handle), intent(in)            :: h
    integer,           intent(out)           :: layout_index
    integer,           intent(out), optional :: ierr

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (h) bind (C, name="esio_field_layout_get")
        import
        integer(c_int)                       :: impl
        type(esio_handle), intent(in), value :: h
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    layout_index = impl(h)
    if (present(ierr)) ierr = 0  ! FIXME: See bug #1178

  end subroutine esio_field_layout_get

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_field_layout_set (h, layout_index, ierr)

    type(esio_handle), intent(in)            :: h
    integer,           intent(in)            :: layout_index
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (h, layout_index) bind (C, name="esio_field_layout_set")
        import
        integer(c_int)                       :: impl
        type(esio_handle), intent(in), value :: h
        integer(c_int),    intent(in), value :: layout_index
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    stat = impl(h, layout_index)
    if (present(ierr)) ierr = stat

  end subroutine esio_field_layout_set

!!@}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_error (reason, file, line, esio_errno)

    character(len=*),  intent(in) :: reason
    character(len=*),  intent(in) :: file
    integer,           intent(in) :: line
    integer,           intent(in) :: esio_errno

#ifndef DOXYGEN_SHOULD_SKIP_THIS
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
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

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
