!#
!#  Program designed to 'chunk' a field into multiple pieces
!#  So we can load it up into an HDF-like storage
!#  Ideally, field should be similar to what is currently 
!#  Used in PSDNS
!#
!#  This routine will therefore be an early test for 
!#  The ESIO library
!#
!#  Lets set the field as follows:
!#  f(x,y,z) = x + y^2 + z^3  
!#
!#  This way, we know exactly how the domain should 'look'
!#  If we were to print it all out. 
!#

program main
  use mpi    ! for mpi

  ! modulefile contains problem size, initial conditions, precision, etc.
  use header, only: ny,nx,nz,nc,xist,xisz,zjst,zjsz,u,initialize_problem,initialize_field
  implicit none

  integer(4) :: myid, numprocs, ierr,a
  integer(4) :: i, leapfrog

  character(len=20) :: filename
  character(len=20) :: dataset
  character(len=20) :: ci
  character(len=20) :: overwrite="o"//char(0)

  real(8) :: stime,etime ! timers

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)
  call mpi_comm_size(mpi_comm_world, numprocs, ierr)

  if(myid.eq.0) write(6,*) "initializing problem"
  call initialize_problem(myid,numprocs) ! subroutine to 'chunk' problem -- header module

  if(myid.eq.0) write(6,*) "initializing field"
  call initialize_field(myid)            ! sets initial field (local field is u)

  ! initial filename for writing field
  filename="field/outpen.1.h5"//char(0)
  dataset="u.field"//char(0) 

  stime=mpi_wtime() !mpi timer
  call esiof_write_field(ny,nx,nz,nc,xist,xisz,zjst,zjsz,u,filename,dataset,overwrite);     
  etime=mpi_wtime() !mpi end timer

  if(myid.eq.0) write(6,*) "program is finalizing, time was: ", etime-stime
  if(myid.eq.0) write(6,*) "In ",leapfrog," iterations"
  if(myid.eq.0) write(6,*) "wrote / read: ", leapfrog*(nc*8*ny*nz*real(nx/2.))/(10e6), " MB"
  if(myid.eq.0) write(6,*) "speed: ", (leapfrog*(nc*8*ny*nz*nx/2.)/(10e6))/(etime-stime), " MB/s"
  call mpi_finalize(ierr)
  
end program main
