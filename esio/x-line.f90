!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! esio 0.2.0: ExaScale IO library for turbulence simulation restart files
!! http://pecos.ices.utexas.edu/
!!
!! Copyright (C) 2010, 2011, 2012 The PECOS Development Team
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
#ifndef VECTORVALUED
  CTYPE,             FINTENT               :: line(*)
#else
  CTYPE,             FINTENT               :: line(1,*)
  integer,           intent(in)            :: ncomponents
#endif
  integer,           intent(in),  optional :: astride
#ifdef HASCOMMENT
  character(len=*),  intent(in),  optional :: comment
#endif
  integer,           intent(out), optional :: ierr

  integer :: stat, tmp_astride
#ifdef HASCOMMENT
  character(len=1,kind=c_char), parameter :: null_char(1) = (/ c_null_char /)
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  interface
    function line_impl (handle, name, line,          &
                        astride                      &
#ifdef VECTORVALUED
                        ,ncomponents                 &
#endif
#ifdef HASCOMMENT
                        ,comment                     &
#endif
                       ) bind (C, name=CBINDNAME)
      import
      integer(c_int)                                  :: line_impl
      type(esio_handle),            intent(in), value :: handle
      character(len=1,kind=c_char), intent(in)        :: name(*)
      CTYPE,                        FINTENT           :: line(*)
      integer(c_int),               intent(in), value :: astride
#ifdef VECTORVALUED
      integer(c_int),               intent(in), value :: ncomponents
#endif
#ifdef HASCOMMENT
      character(len=1,kind=c_char), intent(in)        :: comment(*)
#endif
    end function line_impl
  end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

  tmp_astride = 0
  if (present(astride)) tmp_astride = astride

#ifndef HASCOMMENT
  stat = line_impl(handle, esio_f_c_string(name), line,  &
                   tmp_astride                           &
#ifdef VECTORVALUED
                   ,ncomponents                          &
#endif
                  )
#else  /* HASCOMMENT */
  if (present(comment)) then
    stat = line_impl(handle, esio_f_c_string(name), line,  &
                     tmp_astride                           &
#ifdef VECTORVALUED
                     ,ncomponents                          &
#endif
                     ,esio_f_c_string(comment)             &
                    )
  else
    stat = line_impl(handle, esio_f_c_string(name), line,  &
                     tmp_astride                           &
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
