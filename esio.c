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
 * $Id: esio.c 3743 2009-07-17 18:09:07Z rhys $
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

//#include "config.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include <esio/esio.h>
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
       Might make the most sense to just create a different library for single prec
       and call it 'sesio' and 'desio' for single and double precision

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
  
  note that we are essentially storing 3d dimensional data as a 2-d set
  the data is stored long as y pencils, with nz elements continguous. When the nz 
  wavenumbers are filled, the routine switches to a different nx.

  Thus, the 2d box size is ny x (nz*nx)

  INPUT : string (file name), ny,nx,nz,nc, length of the file passed
  OUTPUT:

************************ */
void esio_init_(char FILE1[], int *ny,int *nx,int *nz,int *nc,int fname_len)
{
  /* this routine needs to be expanded to add another dimension
     where we consider the three components of velocity
  */

  int RANK=2;  /* dimension of storage */
  int DIM1=*ny; // long y-pencils
  int DIM2=(*nx)*(*nz); //rest of the data

  hid_t file1, dataset1;
  hid_t mid1, fid1;
  hsize_t fdim[] = {DIM1, DIM2};
  double buf1[DIM1][DIM2];
  herr_t ret;
  uint  i, j, k;
  
  // for testing -- initialize field
  for ( i = 0; i < DIM1; i++ ) 
    for ( j = 0; j < DIM2; j++ )
      buf1[i][j] = .01;
  
  file1 = H5Fcreate(FILE1, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);                     // create file
  fid1  = H5Screate_simple (RANK, fdim, NULL);                                           // create dataset
  dataset1 = H5Dcreate (file1, "velocity_field", H5T_NATIVE_DOUBLE, fid1, H5P_DEFAULT);  // build dataset
  ret = H5Dwrite (dataset1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf1);     // write field    
  ret = H5Dclose (dataset1);                                                          // close dataset   
  ret = H5Sclose (fid1);                                                              //
  ret = H5Fclose (file1);                                                             // close file     
  
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
          currently just a generic write

  Input: 
  Output:
  
************************ */
void esio_write_()
{
  hid_t file_id,dataset_id;
  herr_t status;
  int i, j; 
  double dset_data[6][4]; // will want to convert this to an array to write
  
  for(int i=0;i<6;i++)
    for(int j=0;j<4;j++)    
      dset_data[i][j]=1;
  
  // open existing file
  file_id = H5Fopen(FILE1, H5F_ACC_RDWR, H5P_DEFAULT);
  
  //open existing dataset for velocity field (should string this out)
  dataset_id = H5Dopen(file_id, "/velocity_field");
  
  // write to proper location in memory (in buffer)
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);
  /* 
     
     The third argument is the memory dataspace.
     The fourth is the file space indentifier
     Specifying H5S_ALL in these slots means that dataspace will be used 
     I doubt we want to do that
  */

  status = H5Dclose(dataset_id); // Close the dataset. 
  status = H5Fclose(file_id);    // Close the file.

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
  hid_t file_id,dataset_id;
  herr_t status;
  double dset_data[6][4]; // will want to convert this to an array to write
  
  for(int i=0;i<6;i++)
    for(int j=0;j<4;j++)    
      dset_data[i][j]=1;
  
  // open existing file
  file_id = H5Fopen(FILE1, H5F_ACC_RDWR, H5P_DEFAULT);
  
  //open existing dataset for velocity field (should string this out)
  dataset_id = H5Dopen(file_id, "/velocity_field");
  
  // write to proper location in memory (in buffer)
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);
  /* 
     
     The third argument is the memory dataspace.
     The fourth is the file space indentifier
     Specifying H5S_ALL in these slots means that dataspace will be used 
     I doubt we want to do that -- prefer to make it a limited dataspace
  */

  status = H5Dclose(dataset_id); // Close the dataset. 
  status = H5Fclose(file_id);    // Close the file.
  
}
/* ************************
  Function esio_read_field

  intent: reads chunk of velocity field to file

  Input: field variable, and location in ny,nx,nz box
  Output:
  
************************ */

void esio_read_field_()
{

  //status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);


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
