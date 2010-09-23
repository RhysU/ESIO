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
  use header ! modulefile contains problem size, initial conditions, precision, etc.
  use esiof  ! this is a module with the ESIO fortran interface
  implicit none

  integer(4) :: myid, numprocs, ierr,a
  integer(4) :: i, leapfrog

  character(len=20) :: filename
  character(len=20) :: ci

  real(8) :: stime,etime ! timers

  leapfrog = 10  ! number of read/write restarts to make

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)
  call mpi_comm_size(mpi_comm_world, numprocs, ierr)

  if(myid.eq.0) write(6,*) "lets get this party started"
  write(6,*) "Task: ", myid, " of ", numprocs," Reporting in"

  if(myid.eq.0) write(6,*) "initializing problem"
  call initialize_problem(myid,numprocs) ! subroutine to 'chunk' problem -- header module

  if(myid.eq.0) write(6,*) "initializing field"
  call initialize_field()                ! sets initial field (local field is u)

  ! initial filename for writing field
  filename="field/outpen.1.h5"//char(0)
  
  stime=mpi_wtime()

  ! allow for multiple writes/reads
  do i=1,leapfrog               

     ! open hdf file -- pass filename
     call esiof_init(filename, ny, nx, nz, nc, myid)
     
     if(myid.eq.0) write(6,*) "write no. ",i
     call esiof_write_field(filename,ny,xist,zjst,xisz,zjsz,nc,u) 
     
     u(:,:,:,:)=0 ! set field to zero to be sure 
     
     if(myid.eq.0) write(6,*) "read no. ",i
     call esiof_read_field(filename,ny,xist,zjst,xisz,zjsz,nc,u)

     write(ci,*) i+1
     ci=adjustl(ci)
     filename="field/outpen."//trim(ci)//".h5"//char(0) ! save new file each timestep
  enddo
  
  etime=mpi_wtime()  

  if(myid.eq.0) write(6,*) "program is finalizing, time was: ", etime-stime, "for ",leapfrog," iterations"
  call mpi_finalize(ierr)
  
end program main
