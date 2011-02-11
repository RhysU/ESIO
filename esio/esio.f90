!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! esio 0.1.2: ExaScale IO library for turbulence simulation restart files
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

! TODO Allow Fortran to detect invalid handle before other failure
! TODO Disable error handling when ierr is present??

!> \file
!! Provides ESIO's Fortran-based public API following the library's
!! \ref concepts "usage concepts".  The Fortran API almost exactly
!! follows ESIO's C-based public API documented within \ref esio.h
!! with only a few exceptions:
!! <ol>
!!   <li>
!!     All Fortran functionality is exposed through subroutines instead
!!     of the C API's functions.
!!   </li>
!!   <li>
!!     The Fortran API uses <tt>integer</tt>, <tt>single</tt>, and
!!     <tt>double</tt> to refer to the types the C API calls <tt>int</tt>,
!!     <tt>float</tt> and <tt>double</tt>.
!!   </li>
!!   <li>
!!     The \c intent(in) and \c intent(out) semantics of each call
!!     are identical to the C version everywhere with the exception
!!     of esio::esio_handle_initialize(), esio::esio_file_path(),
!!     esio::esio_string_get(), and esio::esio_field_layout_get().
!!     These routines have an extra \c intent(out) parameter in Fortran.
!!   </li>
!!   <li>
!!     Logical-like \c int arguments within the C API, for example
!!     esio_file_create()'s \c overwrite and esio_file_open()'s \c readwrite,
!!     have Fortran type \c logical.
!!   </li>
!!   <li>
!!     The line, plane, and field routines place their \c ncomponents,
!!     A, B, and C direction arguments in "Fortran-order"
!!     (that is, fastest index first).
!!   </li>
!!   <li>
!!     All \c astart, \c bstart, and \c cstart arguments are one-indexed
!!     within Fortran as opposed to the C API's zero-based indexing.
!!   </li>
!!   <li>
!!     Stride arguments are optional within Fortran and default to zero.
!!     As in a C API, a zero stride is shorthand for specifying that
!!     the associated data is contiguous in memory.
!!   </li>
!!   <li>
!!     All subroutines contain an optional \c intent(out) \c ierr parameter
!!     that will be zero on successful return and non-zero on failure.
!!     All methods in this module also invoke ESIO's \ref error.h
!!     "error handling mechanisms" on failure.
!!   </li>
!!   <li>
!!     A custom error handler may be registered through the C API in \ref
!!     error.h, but no mechanism is provided to register an error handler
!!     within Fortran.
!!   </li>
!! </ol>

module esio

! Use ISO_C_BINDING for interoperation with ESIO's C implementation.
! Rename C_PTR to ESIO_HANDLE to make the handle object types more opaque.
! State traverses the language boundary as TYPE(ESIO_HANDLE)s.
! Renaming also prevents collision if caller uses ISO_C_BINDING directly.
  use, intrinsic :: iso_c_binding, only: c_char,               &
                                         c_double,             &
                                         c_double_complex,     &
                                         c_float,              &
                                         c_float_complex,      &
                                         c_int,                &
                                         c_null_char,          &
                                         c_associated,         &
                                         esio_handle => c_ptr

  use esio_c_binding, only: esio_f_c_string, esio_f_c_logical, &
                            esio_c_f_stringcopy,               &
                            esio_c_free

  implicit none  ! Nothing implicit

! C interoperation details are kept hidden from the client...
  private :: c_char, c_double, c_double_complex
  private :: c_float, c_float_complex, c_int, c_null_char
  private :: c_associated
  private :: esio_f_c_string, esio_f_c_logical
  private :: esio_c_f_stringcopy
  private :: esio_c_free
! ... with the exception of our opaque handle object type.
  public :: esio_handle

! TODO Allow Fortran to use customizable error handling
! Error handling routine
  private :: esio_error

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Initializing and finalizing a state handle
!! See \ref conceptshandles "handles" for the associated semantics.
!!@{

  subroutine esio_handle_initialize (handle, comm, ierr)

    type(esio_handle), intent(out)           :: handle
    integer,           intent(in)            :: comm
    integer,           intent(out), optional :: ierr

    ! See C routine esio_handle_initialize_fortran
    ! for details on MPI communicator interoperation
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function handle_impl (comm)  &
                           bind (C, name="esio_handle_initialize_fortran")
        import
        type(esio_handle)          :: handle_impl
        integer, intent(in), value :: comm  ! Note integer not integer(c_int)
      end function handle_impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    handle = handle_impl(comm)
    if (present(ierr)) then
      if (c_associated(handle)) then
        ierr = 0  ! 0 == ESIO_SUCCESS
      else
        ierr = 5  ! 5 == ESIO_EFAILED
      end if
    end if

  end subroutine esio_handle_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_handle_finalize (handle, ierr)
    type(esio_handle), intent(in)            :: handle
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function handle_impl (handle)  &
                           bind (C, name="esio_handle_finalize")
        import
        integer(c_int)                       :: handle_impl
        type(esio_handle), intent(in), value :: handle
      end function handle_impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    stat = handle_impl(handle)
    if (present(ierr)) ierr = stat

  end subroutine esio_handle_finalize

!> @}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Opening and closing files
!! See \ref conceptsfiles "file concepts" for more details.
!!@{

  subroutine esio_file_create (handle, file, overwrite, ierr)

    type(esio_handle), intent(in)            :: handle
    character(len=*),  intent(in)            :: file
    logical,           intent(in)            :: overwrite
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function file_impl (handle, file, overwrite)  &
                         bind (C, name="esio_file_create")
        import
        integer(c_int)                                  :: file_impl
        type(esio_handle),            intent(in), value :: handle
        character(len=1,kind=c_char), intent(in)        :: file(*)
        integer(c_int),               intent(in), value :: overwrite
      end function file_impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    stat = file_impl(handle,                       &
                     esio_f_c_string(file),        &
                     esio_f_c_logical(overwrite))
    if (present(ierr)) ierr = stat

  end subroutine esio_file_create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_file_open (handle, file, readwrite, ierr)

    type(esio_handle), intent(in)            :: handle
    character(len=*),  intent(in)            :: file
    logical,           intent(in)            :: readwrite
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function file_impl (handle, file, readwrite)  &
                         bind (C, name="esio_file_open")
        import
        integer(c_int)                                  :: file_impl
        type(esio_handle),            intent(in), value :: handle
        character(len=1,kind=c_char), intent(in)        :: file(*)
        integer(c_int),               intent(in), value :: readwrite
      end function file_impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    stat = file_impl(handle,                       &
                     esio_f_c_string(file),        &
                     esio_f_c_logical(readwrite))
    if (present(ierr)) ierr = stat

  end subroutine esio_file_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_file_clone (handle, srcfile, dstfile, overwrite, ierr)

    type(esio_handle), intent(in)            :: handle
    character(len=*),  intent(in)            :: srcfile
    character(len=*),  intent(in)            :: dstfile
    logical,           intent(in)            :: overwrite
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function file_impl (handle, srcfile, dstfile, overwrite)  &
                         bind (C, name="esio_file_clone")
        import
        integer(c_int)                                  :: file_impl
        type(esio_handle),            intent(in), value :: handle
        character(len=1,kind=c_char), intent(in)        :: srcfile(*)
        character(len=1,kind=c_char), intent(in)        :: dstfile(*)
        integer(c_int),               intent(in), value :: overwrite
      end function file_impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    stat = file_impl(handle,                      &
                     esio_f_c_string(srcfile),    &
                     esio_f_c_string(dstfile),    &
                     esio_f_c_logical(overwrite))
    if (present(ierr)) ierr = stat

  end subroutine esio_file_clone

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_file_path (handle, file_path, ierr)

    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer

    type(esio_handle), intent(in)            :: handle
    character(len=*),  intent(out)           :: file_path
    integer,           intent(out), optional :: ierr

    type(c_ptr)                              :: tmp_p

!   The C implementation returns newly allocated memory
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function file_impl (handle) bind (C, name="esio_file_path")
        import
        type(c_ptr)                          :: file_impl
        type(esio_handle), intent(in), value :: handle
      end function file_impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    tmp_p = file_impl(handle)
    if (esio_c_f_stringcopy(tmp_p, file_path)) then
      if (present(ierr)) ierr = 0  ! 0 == ESIO_SUCCESS
    else
      if (present(ierr)) then
        ierr = 5  ! 5 == ESIO_EFAILED
      else
        call esio_error('esio_file_path failed but ierr was not supplied', &
                        __FILE__, __LINE__, 5)
        call abort
      endif
    end if
    call esio_c_free(tmp_p)

  end subroutine esio_file_path

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_file_flush (handle, ierr)

    type(esio_handle), intent(in)            :: handle
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function file_impl (handle) bind (C, name="esio_file_flush")
        import
        integer(c_int)                       :: file_impl
        type(esio_handle), intent(in), value :: handle
      end function file_impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    stat = file_impl(handle)
    if (present(ierr)) ierr = stat

  end subroutine esio_file_flush

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_file_close (handle, ierr)

    type(esio_handle), intent(in)            :: handle
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function file_impl (handle) bind (C, name="esio_file_close")
        import
        integer(c_int)                       :: file_impl
        type(esio_handle), intent(in), value :: handle
      end function file_impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    stat = file_impl(handle)
    if (present(ierr)) ierr = stat

  end subroutine esio_file_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_file_close_restart (handle, restart_template,  &
                                      retain_count, ierr)

    type(esio_handle), intent(in)            :: handle
    character(len=*),  intent(in)            :: restart_template
    integer,           intent(in)            :: retain_count
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function file_impl (handle, restart_template, retain_count)  &
                         bind (C, name="esio_file_close_restart")
        import
        integer(c_int)                                  :: file_impl
        type(esio_handle),            intent(in), value :: handle
        character(len=1,kind=c_char), intent(in)        :: restart_template(*)
        integer(c_int),               intent(in), value :: retain_count
      end function file_impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    stat = file_impl(handle, esio_f_c_string(restart_template), retain_count)
    if (present(ierr)) ierr = stat

  end subroutine esio_file_close_restart

!!@}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Manipulating string-valued attributes
!! See \ref conceptsattributes "attributes concepts" for more details.
!!@{

  subroutine esio_string_set (handle, location, name, value, ierr)

    type(esio_handle), intent(in)            :: handle
    character(len=*),  intent(in)            :: location
    character(len=*),  intent(in)            :: name
    character(len=*),  intent(in)            :: value
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function string_impl (handle, location, name, value)  &
                           bind (C, name="esio_string_set")
        import
        integer(c_int)                                  :: string_impl
        type(esio_handle),            intent(in), value :: handle
        character(len=1,kind=c_char), intent(in)        :: location(*)
        character(len=1,kind=c_char), intent(in)        :: name(*)
        character(len=1,kind=c_char), intent(in)        :: value(*)
      end function string_impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    stat = string_impl(handle,                    &
                       esio_f_c_string(location), &
                       esio_f_c_string(name),     &
                       esio_f_c_string(value))
    if (present(ierr)) ierr = stat

  end subroutine esio_string_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_string_get (handle, location, name, value, ierr)

    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer

    type(esio_handle), intent(in)            :: handle
    character(len=*),  intent(in)            :: location
    character(len=*),  intent(in)            :: name
    character(len=*),  intent(out)           :: value
    integer,           intent(out), optional :: ierr

    type(c_ptr)                              :: tmp_p

!   The C implementation returns newly allocated memory
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function string_impl (handle, location, name)  &
                           bind (C, name="esio_string_get")
        import
        type(c_ptr)                                     :: string_impl
        type(esio_handle),            intent(in), value :: handle
        character(len=1,kind=c_char), intent(in)        :: location(*)
        character(len=1,kind=c_char), intent(in)        :: name(*)
      end function string_impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    tmp_p = string_impl(handle,                    &
                        esio_f_c_string(location), &
                        esio_f_c_string(name))
    if (esio_c_f_stringcopy(tmp_p, value)) then
      if (present(ierr)) ierr = 0  ! 0 == ESIO_SUCCESS
    else
      if (present(ierr)) then
        ierr = 5  ! 5 == ESIO_EFAILED
      else
        call esio_error('esio_string_get failed but ierr was not supplied', &
                        __FILE__, __LINE__, 5)
        call abort
      endif
    end if
    call esio_c_free(tmp_p)

  end subroutine esio_string_get

!!@}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Manipulating scalar-valued attributes
!! See \ref conceptsattributes "attribute concepts" for more details.
!!@{

subroutine esio_attribute_write_double (handle, location, name, value, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(in)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_attribute_write_double"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_write_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_attribute_write_single (handle, location, name, value, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(in)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_attribute_write_float"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_write_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_attribute_write_integer (handle, location, name, value, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(in)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_attribute_write_int"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_write_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_attribute_read_double (handle, location, name, value, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_attribute_read_double"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_read_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_attribute_read_single (handle, location, name, value, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_attribute_read_float"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_read_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_attribute_read_integer (handle, location, name, value, ierr)
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

subroutine esio_attribute_writev_double (handle, location, name, value, &
                                         ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(in)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_attribute_writev_double"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_writev_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_attribute_writev_single (handle, location, name, value, &
                                         ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(in)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_attribute_writev_float"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_writev_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_attribute_writev_integer (handle, location, name, value, &
                                          ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(in)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_attribute_writev_int"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_writev_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_attribute_readv_double (handle, location, name, value, &
                                        ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_attribute_readv_double"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_readv_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_attribute_readv_single (handle, location, name, value, &
                                        ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_attribute_readv_float"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_readv_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_attribute_readv_integer (handle, location, name, value, &
                                         ncomponents, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_attribute_readv_int"
#  include "attribute.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_attribute_readv_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_attribute_sizev (handle, location, name, ncomponents, ierr)

    type(esio_handle), intent(in)            :: handle
    character(len=*),  intent(in)            :: location
    character(len=*),  intent(in)            :: name
    integer,           intent(out)           :: ncomponents
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    integer(c_int) :: tmp_ncomponents

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (handle, location, name, ncomponents)  &
                    bind (C, name="esio_attribute_sizev")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: handle
        character(len=1,kind=c_char), intent(in)        :: location(*)
        character(len=1,kind=c_char), intent(in)        :: name(*)
        integer(c_int),               intent(inout)     :: ncomponents
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    stat = impl(handle,                    &
                esio_f_c_string(location), &
                esio_f_c_string(name),     &
                tmp_ncomponents)
    ncomponents = tmp_ncomponents
    if (present(ierr)) ierr = stat

  end subroutine esio_attribute_sizev

!!@}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Manipulating decompositions for line, plane, and field operations.
!!@{

  subroutine esio_line_establish (handle, aglobal, astart, alocal, ierr)

    type(esio_handle), intent(in)            :: handle
    integer,           intent(in)            :: aglobal, astart, alocal
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (handle, aglobal, astart, alocal)  &
                    bind (C, name="esio_line_establish")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: handle
        integer(c_int),               intent(in), value :: aglobal, &
                                                           astart,  &
                                                           alocal
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

!   Note conversion from one- to zero-based starting offsets
    stat = impl(handle, aglobal, astart - 1, alocal)
    if (present(ierr)) ierr = stat

  end subroutine esio_line_establish

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_plane_establish (handle, aglobal, astart, alocal,       &
                                           bglobal, bstart, blocal, ierr)

    type(esio_handle), intent(in)            :: handle
    integer,           intent(in)            :: aglobal, astart, alocal
    integer,           intent(in)            :: bglobal, bstart, blocal
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (handle, bglobal, bstart, blocal,     &
                             aglobal, astart, alocal)     &
                    bind (C, name="esio_plane_establish")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: handle
        integer(c_int),               intent(in), value :: bglobal, &
                                                           bstart,  &
                                                           blocal
        integer(c_int),               intent(in), value :: aglobal, &
                                                           astart,  &
                                                           alocal
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

!   Note conversion from one- to zero-based starting offsets
!   Note reordering Fortran's (a, b) to C's (b, a)
    stat = impl(handle, bglobal, bstart - 1, blocal,  &
                        aglobal, astart - 1, alocal)
    if (present(ierr)) ierr = stat

  end subroutine esio_plane_establish

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_field_establish (handle, aglobal, astart, alocal,       &
                                           bglobal, bstart, blocal,       &
                                           cglobal, cstart, clocal, ierr)

    type(esio_handle), intent(in)            :: handle
    integer,           intent(in)            :: aglobal, astart, alocal
    integer,           intent(in)            :: bglobal, bstart, blocal
    integer,           intent(in)            :: cglobal, cstart, clocal
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (handle, cglobal, cstart, clocal,     &
                             bglobal, bstart, blocal,     &
                             aglobal, astart, alocal)     &
                    bind (C, name="esio_field_establish")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: handle
        integer(c_int),               intent(in), value :: cglobal, &
                                                           cstart,  &
                                                           clocal
        integer(c_int),               intent(in), value :: bglobal, &
                                                           bstart,  &
                                                           blocal
        integer(c_int),               intent(in), value :: aglobal, &
                                                           astart,  &
                                                           alocal
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

!   Note conversion from one- to zero-based starting offsets
!   Note reordering Fortran's (a, b, c) to C's (c, b, a)
    stat = impl(handle, cglobal, cstart - 1, clocal,  &
                        bglobal, bstart - 1, blocal,  &
                        aglobal, astart - 1, alocal)
    if (present(ierr)) ierr = stat

  end subroutine esio_field_establish

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_line_established (handle, aglobal, astart, alocal, ierr)

    type(esio_handle), intent(in)            :: handle
    integer,           intent(out)           :: aglobal, astart, alocal
    integer,           intent(out), optional :: ierr

    integer        :: stat
    integer(c_int) :: tmp_aglobal, tmp_astart, tmp_alocal

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (handle, aglobal, astart, alocal)  &
                    bind (C, name="esio_line_established")
        import
        integer(c_int)                                    :: impl
        type(esio_handle),            intent(in),   value :: handle
        integer(c_int),               intent(inout)       :: aglobal, &
                                                             astart,  &
                                                             alocal
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

!   Note conversion from zero- to one-based starting offsets
    stat = impl(handle, tmp_aglobal, tmp_astart, tmp_alocal)
    aglobal = tmp_aglobal
    astart  = tmp_astart + 1
    alocal  = tmp_alocal
    if (present(ierr)) ierr = stat

  end subroutine esio_line_established

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_plane_established (handle, aglobal, astart, alocal,       &
                                             bglobal, bstart, blocal, ierr)

    type(esio_handle), intent(in)            :: handle
    integer,           intent(out)           :: aglobal, astart, alocal
    integer,           intent(out)           :: bglobal, bstart, blocal
    integer,           intent(out), optional :: ierr

    integer        :: stat
    integer(c_int) :: tmp_aglobal, tmp_astart, tmp_alocal
    integer(c_int) :: tmp_bglobal, tmp_bstart, tmp_blocal

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (handle, bglobal, bstart, blocal,  &
                             aglobal, astart, alocal)  &
                    bind (C, name="esio_plane_established")
        import
        integer(c_int)                                    :: impl
        type(esio_handle),            intent(in),   value :: handle
        integer(c_int),               intent(inout)       :: bglobal, &
                                                             bstart,  &
                                                             blocal
        integer(c_int),               intent(inout)       :: aglobal, &
                                                             astart,  &
                                                             alocal
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

!   Note conversion from zero- to one-based starting offsets
!   Note reordering Fortran's (a, b) to C's (b, a)
    stat = impl(handle, tmp_bglobal, tmp_bstart, tmp_blocal,  &
                        tmp_aglobal, tmp_astart, tmp_alocal)
    bglobal = tmp_bglobal
    bstart  = tmp_bstart + 1
    blocal  = tmp_blocal
    aglobal = tmp_aglobal
    astart  = tmp_astart + 1
    alocal  = tmp_alocal
    if (present(ierr)) ierr = stat

  end subroutine esio_plane_established

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_field_established (handle, aglobal, astart, alocal,       &
                                             bglobal, bstart, blocal,       &
                                             cglobal, cstart, clocal, ierr)

    type(esio_handle), intent(in)            :: handle
    integer,           intent(out)           :: aglobal, astart, alocal
    integer,           intent(out)           :: bglobal, bstart, blocal
    integer,           intent(out)           :: cglobal, cstart, clocal
    integer,           intent(out), optional :: ierr

    integer        :: stat
    integer(c_int) :: tmp_aglobal, tmp_astart, tmp_alocal
    integer(c_int) :: tmp_bglobal, tmp_bstart, tmp_blocal
    integer(c_int) :: tmp_cglobal, tmp_cstart, tmp_clocal

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (handle, cglobal, cstart, clocal,  &
                             bglobal, bstart, blocal,  &
                             aglobal, astart, alocal)  &
                    bind (C, name="esio_field_established")
        import
        integer(c_int)                                    :: impl
        type(esio_handle),            intent(in),   value :: handle
        integer(c_int),               intent(inout)       :: cglobal, &
                                                             cstart,  &
                                                             clocal
        integer(c_int),               intent(inout)       :: bglobal, &
                                                             bstart,  &
                                                             blocal
        integer(c_int),               intent(inout)       :: aglobal, &
                                                             astart,  &
                                                             alocal
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

!   Note conversion from zero- to one-based starting offsets
!   Note reordering Fortran's (a, b, c) to C's (c, b, a)
    stat = impl(handle, tmp_cglobal, tmp_cstart, tmp_clocal,  &
                        tmp_bglobal, tmp_bstart, tmp_blocal,  &
                        tmp_aglobal, tmp_astart, tmp_alocal)
    cglobal = tmp_cglobal
    cstart  = tmp_cstart + 1
    clocal  = tmp_clocal
    bglobal = tmp_bglobal
    bstart  = tmp_bstart + 1
    blocal  = tmp_blocal
    aglobal = tmp_aglobal
    astart  = tmp_astart + 1
    alocal  = tmp_alocal
    if (present(ierr)) ierr = stat

  end subroutine esio_field_established

!!@}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Manipulating distributed, one-dimensional, scalar-valued data
!! See \ref conceptslines "line concepts" for more details.
!!@{

subroutine esio_line_write_double (handle, name, line,               &
                                   astride,                          &
                                   comment, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define HASCOMMENT
#  define FINTENT intent(in)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_line_write_double"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_write_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_line_write_single (handle, name, line,               &
                                   astride,                          &
                                   comment, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define HASCOMMENT
#  define FINTENT intent(in)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_line_write_float"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_write_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_line_write_integer (handle, name, line,               &
                                    astride,                          &
                                    comment, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define HASCOMMENT
#  define FINTENT intent(in)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_line_write_int"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_write_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_line_read_double (handle, name, line,               &
                                  astride,                          &
                                  ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_line_read_double"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_read_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_line_read_single (handle, name, line,               &
                                  astride,                          &
                                  ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_line_read_float"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_read_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_line_read_integer (handle, name, line,               &
                                   astride,                          &
                                   ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_line_read_int"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_read_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_line_size (handle, name, aglobal, ierr)

    type(esio_handle), intent(in)            :: handle
    character(len=*),  intent(in)            :: name
    integer,           intent(out)           :: aglobal
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    integer(c_int) :: tmp_a

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (handle, name, aglobal) bind (C, name="esio_line_size")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: handle
        character(len=1,kind=c_char), intent(in)        :: name(*)
        integer(c_int),               intent(inout)     :: aglobal
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    stat = impl(handle, esio_f_c_string(name), tmp_a)
    aglobal = tmp_a
    if (present(ierr)) ierr = stat

  end subroutine esio_line_size

!!@}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Manipulating distributed, one-dimensional, vector-valued data
!! See \ref conceptslines "line concepts" for more details.
!!@{

subroutine esio_line_writev_double (handle, name, line, ncomponents,  &
                                    astride,                          &
                                    comment, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define HASCOMMENT
#  define FINTENT intent(in)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_line_writev_double"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_writev_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_line_writev_single (handle, name, line, ncomponents,  &
                                    astride,                          &
                                    comment, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define HASCOMMENT
#  define FINTENT intent(in)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_line_writev_float"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_writev_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_line_writev_integer (handle, name, line, ncomponents,  &
                                     astride,                          &
                                     comment, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define HASCOMMENT
#  define FINTENT intent(in)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_line_writev_int"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_writev_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_line_readv_double (handle, name, line, ncomponents,  &
                                   astride,                          &
                                   ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_line_readv_double"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_readv_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_line_readv_single (handle, name, line, ncomponents,  &
                                   astride,                          &
                                   ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_line_readv_float"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_readv_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_line_readv_integer (handle, name, line, ncomponents,  &
                                    astride,                          &
                                    ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_line_readv_int"
#  include "line.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_line_readv_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_line_sizev (handle, name, ncomponents, aglobal, ierr)

    type(esio_handle), intent(in)            :: handle
    character(len=*),  intent(in)            :: name
    integer,           intent(out)           :: ncomponents
    integer,           intent(out)           :: aglobal
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    integer(c_int) :: tmp_a, tmp_ncomponents

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (handle, name, aglobal, ncomponents)  &
                    bind (C, name="esio_line_sizev")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: handle
        character(len=1,kind=c_char), intent(in)        :: name(*)
        integer(c_int),               intent(inout)     :: aglobal
        integer(c_int),               intent(inout)     :: ncomponents
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    stat = impl(handle, esio_f_c_string(name), tmp_a, tmp_ncomponents)
    aglobal = tmp_a
    ncomponents = tmp_ncomponents
    if (present(ierr)) ierr = stat

  end subroutine esio_line_sizev

!!@}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Manipulating distributed, two-dimensional, scalar-valued data
!! See \ref conceptsplanes "plane concepts" for more details.
!!@{

subroutine esio_plane_write_double (handle, name, plane,              &
                                    astride, bstride,                 &
                                    comment, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define HASCOMMENT
#  define FINTENT intent(in)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_plane_write_double"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_write_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_plane_write_single (handle, name, plane,              &
                                    astride, bstride,                 &
                                    comment, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define HASCOMMENT
#  define FINTENT intent(in)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_plane_write_float"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_write_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_plane_write_integer (handle, name, plane,              &
                                     astride, bstride,                 &
                                     comment, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define HASCOMMENT
#  define FINTENT intent(in)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_plane_write_int"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_write_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_plane_read_double (handle, name, plane,              &
                                   astride, bstride,                 &
                                   ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_plane_read_double"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_read_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_plane_read_single (handle, name, plane,              &
                                   astride, bstride,                 &
                                   ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_plane_read_float"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_read_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_plane_read_integer (handle, name, plane,              &
                                    astride, bstride,                 &
                                    ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_plane_read_int"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_read_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_plane_size (handle, name, aglobal, bglobal, ierr)

    type(esio_handle), intent(in)            :: handle
    character(len=*),  intent(in)            :: name
    integer,           intent(out)           :: aglobal
    integer,           intent(out)           :: bglobal
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    integer(c_int) :: tmp_b, tmp_a

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (handle, name, bglobal, aglobal)  &
                    bind (C, name="esio_plane_size")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: handle
        character(len=1,kind=c_char), intent(in)        :: name(*)
        integer(c_int),               intent(inout)     :: bglobal
        integer(c_int),               intent(inout)     :: aglobal
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

!   Note reordering Fortran's (a, b) to C's (b, a)
    stat = impl(handle, esio_f_c_string(name), tmp_b, tmp_a)
    aglobal = tmp_a
    bglobal = tmp_b
    if (present(ierr)) ierr = stat

  end subroutine esio_plane_size

!!@}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \name Manipulating distributed, two-dimensional, vector-valued data
!! See \ref conceptsplanes "plane concepts" for more details.
!!@{

subroutine esio_plane_writev_double (handle, name, plane, ncomponents, &
                                     astride, bstride,                 &
                                     comment, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define HASCOMMENT
#  define FINTENT intent(in)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_plane_writev_double"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_writev_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_plane_writev_single (handle, name, plane, ncomponents, &
                                     astride, bstride,                 &
                                     comment, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define HASCOMMENT
#  define FINTENT intent(in)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_plane_writev_float"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_writev_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_plane_writev_integer (handle, name, plane, ncomponents, &
                                      astride, bstride,                 &
                                      comment, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define HASCOMMENT
#  define FINTENT intent(in)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_plane_writev_int"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_writev_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_plane_readv_double (handle, name, plane, ncomponents, &
                                    astride, bstride,                 &
                                    ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_plane_readv_double"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_readv_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_plane_readv_single (handle, name, plane, ncomponents, &
                                    astride, bstride,                 &
                                    ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_plane_readv_float"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_readv_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_plane_readv_integer (handle, name, plane, ncomponents, &
                                     astride, bstride,                 &
                                     ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_plane_readv_int"
#  include "plane.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_plane_readv_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_plane_sizev (handle, name, ncomponents, &
                               aglobal, bglobal, ierr)

    type(esio_handle), intent(in)            :: handle
    character(len=*),  intent(in)            :: name
    integer,           intent(out)           :: ncomponents
    integer,           intent(out)           :: aglobal
    integer,           intent(out)           :: bglobal
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    integer(c_int) :: tmp_b, tmp_a, tmp_ncomponents

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (handle, name, bglobal, aglobal, ncomponents)  &
                    bind (C, name="esio_plane_sizev")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: handle
        character(len=1,kind=c_char), intent(in)        :: name(*)
        integer(c_int),               intent(inout)     :: bglobal
        integer(c_int),               intent(inout)     :: aglobal
        integer(c_int),               intent(inout)     :: ncomponents
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

!   Note reordering Fortran's (a, b) to C's (b, a)
    stat = impl(handle, esio_f_c_string(name), tmp_b, tmp_a, tmp_ncomponents)
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

subroutine esio_field_write_double (handle, name, field,              &
                                    astride, bstride, cstride,        &
                                    comment, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define HASCOMMENT
#  define FINTENT intent(in)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_field_write_double"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_write_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_field_write_single (handle, name, field,              &
                                    astride, bstride, cstride,        &
                                    comment, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define HASCOMMENT
#  define FINTENT intent(in)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_field_write_float"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_write_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_field_write_integer (handle, name, field,              &
                                     astride, bstride, cstride,        &
                                     comment, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define HASCOMMENT
#  define FINTENT intent(in)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_field_write_int"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_write_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_field_read_double (handle, name, field,              &
                                   astride, bstride, cstride,        &
                                   ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_field_read_double"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_read_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_field_read_single (handle, name, field,              &
                                   astride, bstride, cstride,        &
                                   ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_field_read_float"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_read_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_field_read_integer (handle, name, field,              &
                                    astride, bstride, cstride,        &
                                    ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define FINTENT intent(out)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_field_read_int"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_read_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_field_size (handle, name, aglobal, bglobal, cglobal, ierr)

    type(esio_handle), intent(in)            :: handle
    character(len=*),  intent(in)            :: name
    integer,           intent(out)           :: aglobal
    integer,           intent(out)           :: bglobal
    integer,           intent(out)           :: cglobal
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    integer(c_int) :: tmp_c, tmp_b, tmp_a

#ifndef DOXYGEN_SHOULD_SKIP_THIS
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
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

!   Note reordering Fortran's (a, b, c) to C's (c, b, a)
    stat = impl(handle, esio_f_c_string(name), tmp_c, tmp_b, tmp_a)
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

subroutine esio_field_writev_double (handle, name, field, ncomponents, &
                                     astride, bstride, cstride,        &
                                     comment, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define HASCOMMENT
#  define FINTENT intent(in)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_field_writev_double"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_writev_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_field_writev_single (handle, name, field, ncomponents, &
                                     astride, bstride, cstride,        &
                                     comment, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define HASCOMMENT
#  define FINTENT intent(in)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_field_writev_float"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_writev_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_field_writev_integer (handle, name, field, ncomponents, &
                                      astride, bstride, cstride,        &
                                      comment, ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define HASCOMMENT
#  define FINTENT intent(in)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_field_writev_int"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_writev_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_field_readv_double (handle, name, field, ncomponents, &
                                    astride, bstride, cstride,        &
                                    ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE real(c_double)
#  define CBINDNAME "esio_field_readv_double"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_readv_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_field_readv_single (handle, name, field, ncomponents, &
                                    astride, bstride, cstride,        &
                                    ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE real(c_float)
#  define CBINDNAME "esio_field_readv_float"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_readv_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine esio_field_readv_integer (handle, name, field, ncomponents, &
                                     astride, bstride, cstride,        &
                                     ierr)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#  define VECTORVALUED
#  define FINTENT intent(out)
#  define CTYPE integer(c_int)
#  define CBINDNAME "esio_field_readv_int"
#  include "field.f90"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end subroutine esio_field_readv_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_field_sizev (handle, name, ncomponents,             &
                               aglobal, bglobal, cglobal, ierr)

    type(esio_handle), intent(in)            :: handle
    character(len=*),  intent(in)            :: name
    integer,           intent(out)           :: ncomponents
    integer,           intent(out)           :: aglobal
    integer,           intent(out)           :: bglobal
    integer,           intent(out)           :: cglobal
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

    integer(c_int) :: tmp_c, tmp_b, tmp_a, tmp_ncomponents

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (handle, name, cglobal, bglobal, aglobal, ncomponents)  &
                    bind (C, name="esio_field_sizev")
        import
        integer(c_int)                                  :: impl
        type(esio_handle),            intent(in), value :: handle
        character(len=1,kind=c_char), intent(in)        :: name(*)
        integer(c_int),               intent(inout)     :: cglobal
        integer(c_int),               intent(inout)     :: bglobal
        integer(c_int),               intent(inout)     :: aglobal
        integer(c_int),               intent(inout)     :: ncomponents
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

!   Note reordering Fortran's (a, b, c) to C's (c, b, a)
    stat = impl(handle, esio_f_c_string(name),  &
                tmp_c, tmp_b, tmp_a, tmp_ncomponents)
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

  subroutine esio_field_layout_get (handle, layout_index, ierr)

    type(esio_handle), intent(in)            :: handle
    integer,           intent(out)           :: layout_index
    integer,           intent(out), optional :: ierr

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (handle) bind (C, name="esio_field_layout_get")
        import
        integer(c_int)                       :: impl
        type(esio_handle), intent(in), value :: handle
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    layout_index = impl(handle)
    if (present(ierr)) ierr = 0

  end subroutine esio_field_layout_get

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine esio_field_layout_set (handle, layout_index, ierr)

    type(esio_handle), intent(in)            :: handle
    integer,           intent(in)            :: layout_index
    integer,           intent(out), optional :: ierr
    integer                                  :: stat

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    interface
      function impl (handle, layout_index)  &
                    bind (C, name="esio_field_layout_set")
        import
        integer(c_int)                       :: impl
        type(esio_handle), intent(in), value :: handle
        integer(c_int),    intent(in), value :: layout_index
      end function impl
    end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    stat = impl(handle, layout_index)
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

    call c_impl(esio_f_c_string(reason), esio_f_c_string(file), &
                line, esio_errno)

  end subroutine esio_error

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module esio
