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
  public  :: esiof_read, esiof_read_to_var, esiof_write, esiof_timer, esiof_diff, esiof_test
  private ::      ! everything else

!#######################################################
  contains
!#######################################################

    !##############################
    !#
    !# Function esiof_read
    !# 
    !# INPUT : integer
    !# OUTPUT: integer (2x)
    !#
    !##############################
    function esiof_test(testint)
      integer(4), intent(inout) :: testint

      testint = 2*testint

      return testint
    end function esiof_test

    !##############################
    !#
    !# Function esiof_read
    !# 
    !# INPUT : 
    !# OUTPUT:
    !#
    !##############################
    subroutine esiof_read()
      
      call esioc_read
      
    end subroutine esio_read
    
    !##############################
    !#
    !# Function esiof_read_to_var
    !# 
    !# INPUT : 
    !# OUTPUT: 
    !#
    !##############################
    subroutine esiof_read_to_var()
  
      call esioc_read_to_var()

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
    
      call esioc_write()

    end subroutine esiof_write

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
    
      call esioc_timer()

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
    
      call esioc_diff()

    end subroutine esiof_diff

!##########################################    
    
  end module esiof
