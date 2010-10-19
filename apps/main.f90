!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! ESIO - ExaScale IO library
!!
!! ESIO is a restart file library for HPC exascale IO geared especially for
!! turbulence applications.
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

!> Program designed to 'chunk' a field into multiple pieces so we can
!! load it up into an HDF-like storage.  Ideally, field should be similar
!! to what is currently used in PSDNS.
!!
!! This routine will therefore be an early test for the ESIO library.
program main

  use mpi
  use esio

! The header contains problem size, initial conditions, precision, etc.
! The header declares physical space extents but xis{t,z}/zjs{t,z} and
! local storage u are allocated using the wave space extents.
  use header, only: ny, nx, nz, nc,                      & ! Physical space
                    xist, xisz, zjst, zjsz, u,           & ! Wave space
                    initialize_problem, initialize_field
  implicit none

  integer(4)            :: myid, numprocs, ierr, i
  integer(4), parameter :: niter = 1
  type(esio_state)      :: state

  character(len=20) :: filename  = "outpen.1.h5"
  character(len=20) :: fieldname = "u.field"

  real(8) :: stime, etime ! timers

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)
  call mpi_comm_size(mpi_comm_world, numprocs, ierr)

  if (myid.eq.0) write(*,*) "initializing esio"
  state = esio_init(mpi_comm_world)
  ierr  = esio_file_create(state, filename, .TRUE.)

  if (myid.eq.0) write(*,*) "initializing problem"
  call initialize_problem(myid,numprocs)

  if (myid.eq.0) write(*,*) "initializing field:", ny, nz, nx, nc
  call initialize_field(myid)

! Print domain decomposition details on each MPI rank
! do i = 0, numprocs - 1
!   call mpi_barrier(mpi_comm_world, ierr)
!   if (myid.eq.i) then
!     write(*,*) myid, "local field has shape:  ", shape(u)
!     write(*,*) myid, "local xis{tz}, zjs{tz}: ", xist, xisz, zjst, zjsz
!     flush (6)
!   end if
!   call mpi_barrier(mpi_comm_world, ierr)
! end do

  stime = mpi_wtime()
  do i = 1, niter
!   TODO Use generic, precision-agnostic call
!   Zeros in {a,b,c}stride arguments indicate field is contiguous
    ierr = esio_field_write_double(state, fieldname, u, &
                                   ny, 1,    ny,   0,   &
                                   nx, xist, xisz, 0,   &
                                   nz, zjst, zjsz, 0)
  end do
  etime = mpi_wtime()

  if (myid.eq.0) then
    write(*,*) "program is finalizing, time was: ", etime-stime
    write(*,*) "In ", niter, " iterations"
    write(*,*) "wrote / read: ", &
               niter*(nc*8*ny*nz*real(nx/2.))/(10e6), " MB"
    write(*,*) "speed: ", &
               (niter*(nc*8*ny*nz*nx/2.)/(10e6))/(etime-stime), " MB/s"
  end if

  ierr = esio_file_close(state)
  ierr = esio_finalize(state)
  call mpi_finalize(ierr)

end program main
