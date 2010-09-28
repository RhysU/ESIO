! Fortran-specific bindings for ESIO
! These are bindings requiring more than some name mangling

subroutine esiof_write_field(ny,nx,nz,nc,                     &
                             xst,xsz,zst,zsz,                 &
                             data,filename,dataset,overwrite)
  implicit none

  ! inputs -- not to be edited
  integer(4),intent(in) :: ny
  integer(4),intent(in) :: nx
  integer(4),intent(in) :: nz
  integer(4),intent(in) :: nc
  integer(4),intent(in) :: xst
  integer(4),intent(in) :: xsz
  integer(4),intent(in) :: zst
  integer(4),intent(in) :: zsz

  real(8),intent(in),dimension(ny,nx) :: data

  ! inputs -- must be altered to play nice with C
  character(len=20),intent(inout) :: filename
  character(len=20),intent(inout) :: dataset
  character(len=20),intent(inout) :: overwrite

  ! trim to make C happy
  filename = TRIM(filename) //char(0)
  dataset  = TRIM(dataset)  //char(0)
  overwrite= TRIM(overwrite)//char(0)

  ! call the function
  call esiofb_write_double_field(ny,nx,nz,nc,                     &
                                 xst,xsz,zst,zsz,                 &
                                 data,filename,dataset,overwrite)

end subroutine esiof_write_field
