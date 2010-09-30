! Fortran-specific bindings for ESIO
! These are bindings requiring more than some name mangling

subroutine esiof_write_field(ny,nx,nz,nc,                     &
                             xst,xsz,zst,zsz,                 &
                             data,filename,dataset,overwrite)

  use, intrinsic :: iso_c_binding,   only : c_double,    &
                                            c_int,       &
                                            c_null_char

  implicit none

  ! inputs -- not to be edited
  integer(c_int),intent(in) :: ny
  integer(c_int),intent(in) :: nx
  integer(c_int),intent(in) :: nz
  integer(c_int),intent(in) :: nc
  integer(c_int),intent(in) :: xst
  integer(c_int),intent(in) :: xsz
  integer(c_int),intent(in) :: zst
  integer(c_int),intent(in) :: zsz

  real(c_double),intent(in),dimension(ny,nx) :: data

  character(len=20),intent(in) :: filename
  character(len=20),intent(in) :: dataset
  character(len=20),intent(in) :: overwrite

  ! trim to make C happy
  character(len=20) :: filename_c, dataset_c, overwrite_c
  filename_c =TRIM(filename) //c_null_char
  dataset_c  =TRIM(dataset)  //c_null_char
  overwrite_c=TRIM(overwrite)//c_null_char

  ! call the function
  call esiofb_write_double_field(ny,nx,nz,nc,                     &
                                 xst,xsz,zst,zsz,                 &
                                 data,                            &
                                 filename_c,dataset_c,overwrite_c)

end subroutine esiof_write_field
