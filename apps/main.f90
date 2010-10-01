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
  use mpi

  ! The header contains problem size, initial conditions, precision, etc.
  ! The header declares physical space extents but xis{t,z}/zjs{t,z} and
  ! local storage u are allocated using the wave space extents.
  use header, only: ny, nx, nz, nc,                      & ! Physical space
                    xist, xisz, zjst, zjsz, u,           & ! Wave space
                    initialize_problem, initialize_field
  implicit none

  integer(4) :: myid, numprocs, ierr, a
  integer(4) :: i, leapfrog

  character(len=20) :: filename
  character(len=20) :: dataset
  character(len=20) :: ci
  character(len=20) :: overwrite="o"//char(0)

  real(8) :: stime, etime ! timers

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)
  call mpi_comm_size(mpi_comm_world, numprocs, ierr)

  ! subroutine to 'chunk' problem -- header module
  if (myid.eq.0) write(*,*) "initializing problem"
  call initialize_problem(myid,numprocs)

  ! initial field (local field is u)
  if (myid.eq.0) write(*,*) "initializing field"
  call initialize_field(myid)
  if (myid.eq.0) write(*,*) "global field has shape: ", ny, nz, nx, nc

  ! initial filename for writing field
  filename="outpen.1.h5"//char(0)
  dataset="u.field"//char(0)

  stime=mpi_wtime() !mpi timer
  call esiof_write_field(ny, nx, nz, nc, xist, xisz, zjst, zjsz, u, &
                         filename, dataset, overwrite);
  etime=mpi_wtime() !mpi end timer

  ! TODO Leapfrog details broken in output
  if (myid.eq.0) then
    write(*,*) "program is finalizing, time was: ", etime-stime
    write(*,*) "In ",leapfrog," iterations"
    write(*,*) "wrote / read: ", &
               leapfrog*(nc*8*ny*nz*real(nx/2.))/(10e6), " MB"
    write(*,*) "speed: ", &
               (leapfrog*(nc*8*ny*nz*nx/2.)/(10e6))/(etime-stime), " MB/s"
  end if
  call mpi_finalize(ierr)

end program main
