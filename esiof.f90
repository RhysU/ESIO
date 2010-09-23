!--------------------------------------------------------------------------
!
! Copyright (C) 2009 The PSDNS Development Team
!
! Please see https://wiki.ices.utexas.edu/PSDNS for more information.
!
! This file is part of the esio library.
!
! esio is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! esio is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with esio.  If not, see <http://www.gnu.org/licenses/>.
!
!--------------------------------------------------------------------------

!##########################################################################
!#
!#   ESIO-Fortran Library (ESIOF)
!#
!# Fortran module that acts as a front for the .cpp to HDF5 library
!# User should only have to add line:
!# 'use esio' 
!# to grab access to esio subroutines and functions
!#
!# This module should also layout the fortran API for ESIO
!# As nearly all the function calls are just passes to C routines
!#
!##########################################################################
module esiof
  implicit none

  ! global module variables here

  ! permissions
  public  :: esiof_read, esiof_read_to_var, esiof_write, esiof_timer, esiof_diff
  public  :: esiof_test, esiof_init, esiof_write_field, esiof_read_field

  integer(4) :: precision  !precision in use

!#######################################################
  contains
!#######################################################

    !##############################
    !#
    !# Function esiof_init
    !# 
    !# INPUT : metadata,filename, problem size
    !# OUTPUT: 
    !#
    !##############################
    subroutine esiof_init(filename,ny,nx,nz,nc,np)
      integer(4), intent(in)        :: ny,nx,nz,nc,np
      character(len=20),intent(out) :: filename     

      call esio_init(filename,ny,nx,nz,nc,np)

    end subroutine esiof_init

    !##############################
    !#
    !# Function esiof_read
    !# 
    !# INPUT : integer
    !# OUTPUT: integer (2x)
    !#
    !##############################
    subroutine esiof_test(testint)
      integer(4), intent(inout) :: testint
            
      call esio_test(testint)
      
    end subroutine esiof_test

    !##############################
    !#
    !# Function esiof_read
    !# 
    !# INPUT : 
    !# OUTPUT:
    !#
    !##############################
    subroutine esiof_read()

      call esio_read
      
    end subroutine esiof_read
    
    !##############################
    !#
    !# Function esiof_read_to_var
    !# 
    !# INPUT : 
    !# OUTPUT: 
    !#
    !##############################
    subroutine esiof_read_to_var()
  
      call esio_read_to_var()

    end subroutine esiof_read_to_var

    !##############################
    !#
    !# Function esiof_write
    !# 
    !# INPUT : 
    !# OUTPUT: 
    !#
    !##############################
    subroutine esiof_write()
    
      call esio_write()

    end subroutine esiof_write

    !##############################
    !#
    !# Function esiof_write_field
    !# 
    !# Does not need to know what processor
    !# Just where the data goes (xist,zjst,xisz...)
    !#
    !# INPUT : 
    !# OUTPUT: 
    !#
    !##############################
    subroutine esiof_write_field(filename,ny,xist,zjst,xisz,zjsz,nc,upencil)
      integer(4), intent(in)                           :: ny,xist,zjst,xisz,zjsz,nc
      character(len=20),intent(in)                     :: filename
      complex(8),dimension(ny,zjst,xist,nc),intent(in) :: upencil

      !this is where we need to switch between different schemes-- i.e. hdf, netcdf, etc.
      call esio_write_field(filename,ny,xist,zjst,xisz,zjsz,nc)!,upencil)

    end subroutine esiof_write_field

    !##############################
    !#
    !# Function esiof_read_field
    !# 
    !# INPUT : 
    !# OUTPUT: 
    !#
    !##############################
    subroutine esiof_read_field(filename,ny,xist,zjst,xisz,zjsz,nc,upencil)
      integer(4), intent(in)                            :: ny,xist,zjst,xisz,zjsz,nc
      character(len=20),intent(in)                      :: filename
      complex(8),dimension(ny,zjst,xist,nc)             :: upencil 
!      complex(8),dimension(ny,zjst,xist,nc),intent(out) :: upencil 
            
      call esio_read_field(filename,ny,xist,zjst,xisz,zjsz,nc)!,upencil)

    end subroutine esiof_read_field

    !##############################
    !#
    !# Function esiof_timer
    !# 
    !# When called, return time
    !# Acts in a manner that provides
    !# timing information
    !#
    !# INPUT : 
    !# OUTPUT: 
    !#
    !##############################
    subroutine esiof_timer()

      call esio_timer()

    end subroutine esiof_timer    

    !##############################
    !#
    !# Function esiof_diff
    !#
    !# Finds the difference between 
    !# two esio files
    !# 
    !# INPUT : 
    !# OUTPUT: 
    !#
    !##############################
    subroutine esiof_diff()

      call esio_diff()

    end subroutine esiof_diff

!##########################################
    
  end module esiof
