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

#if !(defined(FINTENT) && defined(CTYPE) && defined(CBINDNAME) && defined(IMPL))
#error "One of FINTENT, CTYPE, CBINDNAME, or IMPL not defined"
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
    function IMPL (handle, name, line,          &
                   astride                      &
#ifdef VECTORVALUED
                   ,ncomponents                 &
#endif
#ifdef HASCOMMENT
                   ,comment                     &
#endif
                  ) bind (C, name=CBINDNAME)
      import
      integer(c_int)                                  :: IMPL
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
    end function IMPL
  end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

  tmp_astride = 0
  if (present(astride)) tmp_astride = astride

#ifndef HASCOMMENT
  stat = IMPL(handle, esio_f_c_string(name), line,  &
              tmp_astride                           &
#ifdef VECTORVALUED
              ,ncomponents                          &
#endif
             )
#else  /* HASCOMMENT */
  if (present(comment)) then
    stat = IMPL(handle, esio_f_c_string(name), line,  &
                tmp_astride                           &
#ifdef VECTORVALUED
                ,ncomponents                          &
#endif
                ,esio_f_c_string(comment)             &
               )
  else
    stat = IMPL(handle, esio_f_c_string(name), line,  &
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
#undef IMPL
#ifdef HASCOMMENT
#undef HASCOMMENT
#endif
#ifdef VECTORVALUED
#undef VECTORVALUED
#endif
