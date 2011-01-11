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

! Designed to be #included from esio.f90 within a subroutine declaration

#if !defined(FINTENT) || !defined(CTYPE) || !defined(CBINDNAME)
#error "One of FINTENT, CTYPE, or CBINDNAME not defined"
#endif

  type(esio_handle), intent(in)            :: handle
  character(len=*),  intent(in)            :: name
  CTYPE,             FINTENT               :: field(*)
#ifdef VECTORVALUED
  integer,           intent(in)            :: ncomponents
#endif
  integer,           intent(in),  optional :: astride
  integer,           intent(in),  optional :: bstride
  integer,           intent(in),  optional :: cstride
  integer,           intent(out), optional :: ierr

  integer :: stat, tmp_astride, tmp_bstride, tmp_cstride

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  interface
    function field_impl (handle, name, field,        &
                         cstride, bstride, astride   &
#ifdef VECTORVALUED
                         ,ncomponents                &
#endif
                        ) bind (C, name=CBINDNAME)
      import
      integer(c_int)                                  :: field_impl
      type(esio_handle),            intent(in), value :: handle
      character(len=1,kind=c_char), intent(in)        :: name(*)
      CTYPE,                        FINTENT           :: field(*)
      integer(c_int),               intent(in), value :: cstride
      integer(c_int),               intent(in), value :: bstride
      integer(c_int),               intent(in), value :: astride
#ifdef VECTORVALUED
      integer(c_int),               intent(in), value :: ncomponents
#endif
    end function field_impl
  end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

  tmp_astride = 0
  if (present(astride)) tmp_astride = astride
  tmp_bstride = 0
  if (present(bstride)) tmp_bstride = bstride
  tmp_cstride = 0
  if (present(cstride)) tmp_cstride = cstride

! Note reordering Fortran's (a, b, c) to C's (c, b, a)
  stat = field_impl(handle, esio_f_c_string(name), field, &
                    tmp_cstride, tmp_bstride, tmp_astride &
#ifdef VECTORVALUED
                    ,ncomponents                          &
#endif
                   )
  if (present(ierr)) ierr = stat

#undef FINTENT
#undef CTYPE
#undef CBINDNAME
#ifdef VECTORVALUED
#undef VECTORVALUED
#endif
