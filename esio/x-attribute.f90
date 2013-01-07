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
!! by the Free Software Foundation, either version 2.1 of the License, or
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

! Designed to be #included from esio.f90 within a subroutine declaration

#if !defined(FINTENT) || !defined(CTYPE) || !defined(CBINDNAME)
#error "One of FINTENT, CTYPE, or CBINDNAME not defined"
#endif

  type(esio_handle), intent(in) :: handle
  character(len=*),  intent(in) :: location
  character(len=*),  intent(in) :: name
#ifndef VECTORVALUED
  CTYPE,             FINTENT    :: value
#else
  CTYPE,             FINTENT    :: value(*)
  integer,           intent(in) :: ncomponents
#endif
  integer,           intent(out), optional :: ierr
  integer                                  :: stat

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  interface
    function attribute_impl (handle, location, name, value    &
#ifdef VECTORVALUED
                             ,ncomponents                     &
#endif
                            ) bind (C, name=CBINDNAME)
      import
      integer(c_int)                                  :: attribute_impl
      type(esio_handle),            intent(in), value :: handle
      character(len=1,kind=c_char), intent(in)        :: location(*)
      character(len=1,kind=c_char), intent(in)        :: name(*)
#ifndef VECTORVALUED
      CTYPE,                        FINTENT           :: value
#else
      CTYPE,                        FINTENT           :: value(*)
      integer(c_int),               intent(in), value :: ncomponents
#endif
    end function attribute_impl
  end interface
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

  stat = attribute_impl(handle,                               &
                        esio_f_c_string(location),            &
                        esio_f_c_string(name),                &
                        value  &
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
