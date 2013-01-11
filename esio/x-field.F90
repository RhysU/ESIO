!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! ESIO 0.1.8: ExaScale IO library for turbulence simulation restart files
!! http://red.ices.utexas.edu/projects/esio/
!!
!! Copyright (C) 2010, 2011, 2012, 2013 The PECOS Development Team
!!
!! This file is part of ESIO.
!!
!! ESIO is free software: you can redistribute it and/or modify
!! it under the terms of the GNU Lesser General Public License as published
!! by the Free Software Foundation, either version 3.0 of the License, or
!! (at your option) any later version.
!!
!! ESIO is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public License
!! along with ESIO.  If not, see <http://www.gnu.org/licenses/>.
!!
!!-----------------------------------------------------------------------el-
!! $Id$

! Designed to be #included from esio.F90 within a subroutine declaration

#if !defined(FINTENT) || !defined(CTYPE) || !defined(CBINDNAME)
#error "One of FINTENT, CTYPE, or CBINDNAME not defined"
#endif

! Permit specifying interface names to workaround, e.g. XLF, bugs.
#ifndef FIELD_IMPL
#define FIELD_IMPL field_impl
#endif

  type(esio_handle), intent(in)            :: handle
  character(len=*),  intent(in)            :: name
#ifndef VECTORVALUED
  CTYPE,             FINTENT               :: field(1,1,*)
#else
  CTYPE,             FINTENT               :: field(1,1,1,*)
  integer,           intent(in)            :: ncomponents
#endif
  integer,           intent(in),  optional :: astride
  integer,           intent(in),  optional :: bstride
  integer,           intent(in),  optional :: cstride
#ifdef HASCOMMENT
  character(len=*),  intent(in),  optional :: comment
#endif
  integer,           intent(out), optional :: ierr

  integer :: stat, tmp_astride, tmp_bstride, tmp_cstride
#ifdef HASCOMMENT
  character(len=1,kind=c_char), parameter :: null_char(1) = (/ c_null_char /)
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  interface
    function FIELD_IMPL (handle, name, field,        &
                         cstride, bstride, astride   &
#ifdef VECTORVALUED
                         ,ncomponents                &
#endif
#ifdef HASCOMMENT
                        ,comment                     &
#endif
                        ) bind (C, name=CBINDNAME)
      import
      integer(c_int)                                  :: FIELD_IMPL
      type(esio_handle),            intent(in), value :: handle
      character(len=1,kind=c_char), intent(in)        :: name(*)
      CTYPE,                        FINTENT           :: field(*)
      integer(c_int),               intent(in), value :: cstride
      integer(c_int),               intent(in), value :: bstride
      integer(c_int),               intent(in), value :: astride
#ifdef VECTORVALUED
      integer(c_int),               intent(in), value :: ncomponents
#endif
#ifdef HASCOMMENT
      character(len=1,kind=c_char), intent(in)        :: comment(*)
#endif
    end function FIELD_IMPL
  end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

  tmp_astride = 0
  if (present(astride)) tmp_astride = astride
  tmp_bstride = 0
  if (present(bstride)) tmp_bstride = bstride
  tmp_cstride = 0
  if (present(cstride)) tmp_cstride = cstride

! Note reordering Fortran's (a, b, c) to C's (c, b, a)
#ifndef HASCOMMENT
  stat = FIELD_IMPL(handle, esio_f_c_string(name), field, &
                    tmp_cstride, tmp_bstride, tmp_astride &
#ifdef VECTORVALUED
                    ,ncomponents                          &
#endif
                   )
#else /* HASCOMMENT */
  if (present(comment)) then
    stat = FIELD_IMPL(handle, esio_f_c_string(name), field, &
                      tmp_cstride, tmp_bstride, tmp_astride &
#ifdef VECTORVALUED
                      ,ncomponents                          &
#endif
                      ,esio_f_c_string(comment)             &
                     )
  else
    stat = FIELD_IMPL(handle, esio_f_c_string(name), field, &
                      tmp_cstride, tmp_bstride, tmp_astride &
#ifdef VECTORVALUED
                      ,ncomponents                          &
#endif
                      ,null_char                            &
                     )
  end if
#endif /* HASCOMMENT */

  if (present(ierr)) ierr = stat

#undef FINTENT
#undef CTYPE
#undef CBINDNAME
#ifdef HASCOMMENT
#undef HASCOMMENT
#endif
#ifdef VECTORVALUED
#undef VECTORVALUED
#endif
#undef FIELD_IMPL
