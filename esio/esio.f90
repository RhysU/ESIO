!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! esio 0.0.1: ExaScale IO library for turbulence simulation restart files
!! 
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

! TODO Document Fortran API
! TODO Allow Fortran to detect invalid state before other failure
! TODO Allow Fortran to use customizable error handling

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
                                         esio_state => c_ptr

  implicit none  ! Nothing implicit
  private        ! Everything default private

! C interoperation details are kept hidden from the client...
! ... with the exception of our opaque handle object type.
  public :: esio_state

! Public API
  public :: esio_init, esio_finalize
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

  type(esio_state) function esio_init (comm)

    integer, intent(in) :: comm

    ! See C routine esio_init_fortran re: MPI communicator interoperation
    interface
      function impl (comm) bind (C, name="esio_init_fortran")
        import
        type(esio_state)           :: impl
        integer, intent(in), value :: comm  ! Note integer not integer(c_int)
      end function impl
    end interface

    esio_init = impl(comm)

  end function esio_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_finalize (state)

    type(esio_state), intent(in) :: state

    interface
      function impl (state) bind (C, name="esio_finalize")
        import
        integer(c_int)                      :: impl
        type(esio_state), intent(in), value :: state
      end function impl
    end interface

    esio_finalize = impl(state)

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

  integer function esio_layout_get (state)

    type(esio_state), intent(in) :: state

    interface
      function impl (state) bind (C, name="esio_layout_get")
        import
        integer(c_int)                      :: impl
        type(esio_state), intent(in), value :: state
      end function impl
    end interface

    esio_layout_get = impl(state)

  end function esio_layout_get

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_layout_set (state, layout_tag)

    type(esio_state), intent(in) :: state
    integer,          intent(in) :: layout_tag

    interface
      function impl (state, layout_tag) bind (C, name="esio_layout_set")
        import
        integer(c_int)                      :: impl
        type(esio_state), intent(in), value :: state
        integer(c_int),   intent(in), value :: layout_tag
      end function impl
    end interface

    esio_layout_set = impl(state, layout_tag)

  end function esio_layout_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_file_create (state, file, overwrite)

    type(esio_state), intent(in) :: state
    character(len=*), intent(in) :: file
    logical,          intent(in) :: overwrite

    interface
      function impl (state, file, overwrite) bind (C, name="esio_file_create")
        import
        integer(c_int)                                  :: impl
        type(esio_state),             intent(in), value :: state
        character(len=1,kind=c_char), intent(in)        :: file(*)
        integer(c_int),               intent(in), value :: overwrite
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
      function impl (state, file, readwrite) bind (C, name="esio_file_open")
        import
        integer(c_int)                                  :: impl
        type(esio_state),             intent(in), value :: state
        character(len=1,kind=c_char), intent(in)        :: file(*)
        integer(c_int),               intent(in), value :: readwrite
      end function impl
    end interface

    esio_file_open = impl(state, f_c_string(file), f_c_logical(readwrite))

  end function esio_file_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_file_close (state)

    type(esio_state), intent(in) :: state

    interface
      function impl (state) bind (C, name="esio_file_close")
        import
        integer(c_int)                      :: impl
        type(esio_state), intent(in), value :: state
      end function impl
    end interface

    esio_file_close = impl(state)

  end function esio_file_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_field_size (state, name, aglobal, bglobal, cglobal)

    type(esio_state), intent(in)  :: state
    character(len=*), intent(in)  :: name
    integer,          intent(out) :: aglobal
    integer,          intent(out) :: bglobal
    integer,          intent(out) :: cglobal

    integer(c_int) :: tmp_c, tmp_b, tmp_a

    interface
      function impl (state, name, cglobal, bglobal, aglobal)  &
                    bind (C, name="esio_field_size")
        import
        integer(c_int)                                  :: impl
        type(esio_state),             intent(in), value :: state
        character(len=1,kind=c_char), intent(in)        :: name(*)
        integer(c_int),               intent(inout)     :: cglobal
        integer(c_int),               intent(inout)     :: bglobal
        integer(c_int),               intent(inout)     :: aglobal
      end function impl
    end interface

!   Note reordering Fortran's (a, b, c) to C's (c, b, a)
    esio_field_size = impl(state, f_c_string(name), tmp_c, tmp_b, tmp_a)
    aglobal = tmp_a
    bglobal = tmp_b
    cglobal = tmp_c

  end function esio_field_size


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_field_write_double (state, name, field,               &
                                            aglobal, astart, alocal, astride, &
                                            bglobal, bstart, blocal, bstride, &
                                            cglobal, cstart, clocal, cstride)

    type(esio_state), intent(in) :: state
    character(len=*), intent(in) :: name
    real(c_double),   intent(in) :: field(*)
    integer,          intent(in) :: aglobal, astart, alocal, astride
    integer,          intent(in) :: bglobal, bstart, blocal, bstride
    integer,          intent(in) :: cglobal, cstart, clocal, cstride

    interface
      function impl (state, name, field,               &
                     cglobal, cstart, clocal, cstride, &
                     bglobal, bstart, blocal, bstride, &
                     aglobal, astart, alocal, astride) &
                     bind (C, name="esio_field_write_double")
        import
        integer(c_int)                                  :: impl
        type(esio_state),             intent(in), value :: state
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
    esio_field_write_double = impl(state, f_c_string(name), field,       &
                                   cglobal, cstart - 1, clocal, cstride, &
                                   bglobal, bstart - 1, blocal, bstride, &
                                   aglobal, astart - 1, alocal, astride)

  end function esio_field_write_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_field_write_single (state, name, field,               &
                                            aglobal, astart, alocal, astride, &
                                            bglobal, bstart, blocal, bstride, &
                                            cglobal, cstart, clocal, cstride)

    type(esio_state), intent(in) :: state
    character(len=*), intent(in) :: name
    real(c_float),    intent(in) :: field(*)
    integer,          intent(in) :: aglobal, astart, alocal, astride
    integer,          intent(in) :: bglobal, bstart, blocal, bstride
    integer,          intent(in) :: cglobal, cstart, clocal, cstride

    interface
      function impl (state, name, field,               &
                     cglobal, cstart, clocal, cstride, &
                     bglobal, bstart, blocal, bstride, &
                     aglobal, astart, alocal, astride) &
                     bind (C, name="esio_field_write_float")
        import
        integer(c_int)                                  :: impl
        type(esio_state),             intent(in), value :: state
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
    esio_field_write_single = impl(state, f_c_string(name), field,       &
                                   cglobal, cstart - 1, clocal, cstride, &
                                   bglobal, bstart - 1, blocal, bstride, &
                                   aglobal, astart - 1, alocal, astride)

  end function esio_field_write_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_field_read_double (state, name, field,               &
                                           aglobal, astart, alocal, astride, &
                                           bglobal, bstart, blocal, bstride, &
                                           cglobal, cstart, clocal, cstride)

    type(esio_state), intent(in)  :: state
    character(len=*), intent(in)  :: name
    real(c_double),   intent(out) :: field(*)
    integer,          intent(in)  :: aglobal, astart, alocal, astride
    integer,          intent(in)  :: bglobal, bstart, blocal, bstride
    integer,          intent(in)  :: cglobal, cstart, clocal, cstride

    interface
      function impl (state, name, field,               &
                     cglobal, cstart, clocal, cstride, &
                     bglobal, bstart, blocal, bstride, &
                     aglobal, astart, alocal, astride) &
                     bind (C, name="esio_field_read_double")
        import
        integer(c_int)                                  :: impl
        type(esio_state),             intent(in), value :: state
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
    esio_field_read_double = impl(state, f_c_string(name), field,       &
                                  cglobal, cstart - 1, clocal, cstride, &
                                  bglobal, bstart - 1, blocal, bstride, &
                                  aglobal, astart - 1, alocal, astride)

  end function esio_field_read_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function esio_field_read_single (state, name, field,                &
                                           aglobal, astart, alocal, astride,  &
                                           bglobal, bstart, blocal, bstride,  &
                                           cglobal, cstart, clocal, cstride)

    type(esio_state), intent(in)  :: state
    character(len=*), intent(in)  :: name
    real(c_float),    intent(out) :: field(*)
    integer,          intent(in)  :: aglobal, astart, alocal, astride
    integer,          intent(in)  :: bglobal, bstart, blocal, bstride
    integer,          intent(in)  :: cglobal, cstart, clocal, cstride

    interface
      function impl (state, name, field,               &
                     cglobal, cstart, clocal, cstride, &
                     bglobal, bstart, blocal, bstride, &
                     aglobal, astart, alocal, astride) &
                     bind (C, name="esio_field_read_float")
        import
        integer(c_int)                                  :: impl
        type(esio_state),             intent(in), value :: state
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
    esio_field_read_single = impl(state, f_c_string(name), field,       &
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
