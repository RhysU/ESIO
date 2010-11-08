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

! Designed to be #included from esio.f90
! If VECTORVALUED has been #defined, then ncomponents is added to arguments
#if !defined(FNAME) || !defined(FINTENT) || !defined(CTYPE) || !defined(CBINDNAME)
#error "One of FNAME, FINTENT, CTYPE, or CBINDNAME not defined"
#endif

subroutine FNAME (h, name, plane,                   &
                  aglobal, astart, alocal, astride, &
                  bglobal, bstart, blocal, bstride, &
#ifdef VECTORVALUED
                  ncomponents,                      &
#endif
                  ierr)

  type(esio_handle), intent(in) :: h
  character(len=*),  intent(in) :: name
  CTYPE,             FINTENT    :: plane(*)
  integer,           intent(in) :: aglobal, astart, alocal, astride
  integer,           intent(in) :: bglobal, bstart, blocal, bstride
#ifdef VECTORVALUED
  integer,           intent(in) :: ncomponents
#endif
  integer,           intent(out), optional :: ierr
  integer                                  :: stat

  interface
    function impl (h, name, plane,                   &
                   bglobal, bstart, blocal, bstride, &
                   aglobal, astart, alocal, astride  &
#ifdef VECTORVALUED
                   ,ncomponents                      &
#endif
                   ) bind (C, name=CBINDNAME)
      import
      integer(c_int)                                  :: impl
      type(esio_handle),            intent(in), value :: h
      character(len=1,kind=c_char), intent(in)        :: name(*)
      CTYPE,                        FINTENT           :: plane(*)
      integer(c_int),               intent(in), value :: bglobal, &
                                                         bstart,  &
                                                         blocal,  &
                                                         bstride
      integer(c_int),               intent(in), value :: aglobal, &
                                                         astart,  &
                                                         alocal,  &
                                                         astride
#ifdef VECTORVALUED
      integer(c_int),               intent(in), value :: ncomponents
#endif
    end function impl
  end interface

! Note conversion from one- to zero-based starting offsets
! Note reordering Fortran's (a, b) to C's (b, a)
  stat = impl(h, f_c_string(name), plane,           &
              bglobal, bstart - 1, blocal, bstride, &
              aglobal, astart - 1, alocal, astride  &
#ifdef VECTORVALUED
              ,ncomponents                          &
#endif
             )
  if (present(ierr)) ierr = stat

end subroutine FNAME

#undef FNAME
#undef FINTENT
#undef CTYPE
#undef CBINDNAME
#ifdef VECTORVALUED
#undef VECTORVALUED
#endif
