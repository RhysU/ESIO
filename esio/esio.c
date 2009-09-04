/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PSDNS Development Team
 *
 * Please see https://wiki.ices.utexas.edu/PSDNS for more information.
 *
 * This file is part of the esio library.
 *
 * esio is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * esio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with esio.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * esio.c: esio public api implementations
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include "config.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <esio/esio.h>
#include <hdf5.h>

/* 
   TODO:
    - expand C API to more functions
    - fix up fortran interface
    - build the HDF interface

    Notes:
    - Need cpp guru to look over, check this is an efficient structure

    - I considered making this all into a class, but then we would need
      an object of that class before we call the methods

    - Instead, just went with a whole bunch of functions/subroutines to be compiled
    
    - Remember that Fortran passess all variable by reference, not by value.
      Thus, all incoming variables should be assumed to be pointers
         If we need to override this, we should just use %val() in the fortran code
	 Which forces a pass by value
	 
    Problems:
    - we need to be able to handle double or single precision


 */

/* ************************
  Function esio_test
  
  intent: test function for debugging

  INPUT : integer
  OUTPUT: integer (2x)

************************ */
void esio_test_(int *passed)
{
  printf("\nthis is what I got: %i\n",*passed);
  *passed = 2*(*passed);
}

/* ************************
  Function esio_init
  
  intent: open file, put in metadata, etc.

  INPUT : string (file name) 
  OUTPUT:

************************ */
void esio_init_(char FILE1[], int *ny,int *nx,int *nz,int *nc,int *passed,int fname_len)
{
  int RANK=1;
  int DIM1=*ny;
  int DIM2=*nx;
  int DIM3=*nz;

  hid_t file1, dataset1;
  hid_t mid1, fid1;
  hsize_t fdim[] = {DIM1, DIM2, DIM3};
  int buf1[DIM1][DIM2][DIM3];
  herr_t ret;
  unsigned int  i, j, k;
  
  /* for testing */
  for ( i = 0; i < DIM1; i++ ) 
    for ( j = 0; j < DIM2; j++ )
      for ( k = 0; k < DIM3; k++ )
	buf1[i][j][k] = *passed;
  
  file1 = H5Fcreate(FILE1, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);     
  fid1  = H5Screate_simple (RANK, fdim, NULL);
  dataset1 = H5Dcreate (file1, "Velocity Field", H5T_NATIVE_INT, fid1, H5P_DEFAULT);  /* build dataset */
  ret = H5Dwrite (dataset1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf1);     /* write field   */
  ret = H5Dclose (dataset1);                                                          /* close dataset */
  ret = H5Sclose (fid1);                                                              
  ret = H5Fclose (file1);                                                             /* close file    */
  
}

/* ************************
  Function esio_read

  intent: reads a location

  Input: 
  Output:
  
************************ */
void esio_read_()
{

}


/* ************************
  Function esio_read_to_var
  
  intent: populates a given variable from file

  Input: 
  Output:
  
************************ */
void esio_read_to_var_()
{

}


/* ************************
  Function esio_write

  intent: writes a variable to file

  Input: 
  Output:
  
************************ */
void esio_write_()
{
  printf("\nIn the c routine esio_write, printing\n");
}


/* ****************************************
  Function esio_write_field

  intent: writes chunk of velocity field to file

  Input: field variable, and location in ny,nx,nz box
  Output:
  
  **************************************** */

void esio_write_field_(char FILE1[], int *ny,int *xist,int *zjst, int *xisz, int *zjsz, int fname_len)
{


}
/* ************************
  Function esio_read_field

  intent: reads chunk of velocity field to file

  Input: field variable, and location in ny,nx,nz box
  Output:
  
************************ */

void esio_read_field_()
{

}


/* ************************
  Function esio_timer

  intent: provides timing information

  Input: 
  Output:
  
************************ */
void esio_timer_()
{

}


/* ************************
  Function esio_diff

  Input: 
  Output:
  
************************ */

void esio_diff_()
{

}
