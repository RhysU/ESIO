!!-----------------------------------------------------------------------bl-
!!--------------------------------------------------------------------------
!!
!! esio 0.1.1: ExaScale IO library for turbulence simulation restart files
!! http://pecos.ices.utexas.edu/
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
!! $Id$

!> Module that contains the problem in use.  Parallel problem
!! decomposition based upon P3DFFT's logic but may differ in
!! subtle ways.
!!
!! The problem defines a box with lengths ny,nx,nz
!! with indices over 1 -> ny+1, etc.
!!
!! The variables one likely needs:
!! u (field), xist,xien,zjst,zjen (coordinates of field)
module p3dfft_like_header

  use mpi
  implicit none

  integer(4), parameter :: nx=32                 ! physical space lengths
  integer(4), parameter :: ny=1000
  integer(4), parameter :: nz=16
  integer(4), parameter :: nc=1                  ! number of scalar components

  integer(4) :: xist,xjst,xisz,xjsz,xien,xjen    ! pencil widths and heights
  integer(4) :: yist,yjst,yisz,yjsz,yien,yjen
  integer(4) :: zist,zjst,zisz,zjsz,zien,zjen

  real(8),dimension(:,:,:,:),allocatable :: u    ! field: u(ny,zjsz,xisz,nc)

! Public declarations
  public :: u, xisz, zjsz, ny, nc                ! velocity field data - local
  public :: xist,xien,zjst,zjen                  ! absolute field coordinates
  public :: nz,nx                                ! field size
  public :: initialize_problem, initialize_field ! public subroutines

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> This routine generates a local initial field for each MPI rank.
!!
!! For reference, xisz (x-size), zjsz, and ny are local coordinates
!! on each processor.  Global sizes are ny, nx, and nz.  Absolute
!! coordinates are given by xist xien, zjst zjen, and 1 ny in the
!! x, z, and y directions.
  subroutine initialize_field(proc)

    integer(4) :: i,j,k,n,proc

    !this loop is stride-1: u(ny,zjsz,xisz)
    do n=1,nc
      do i=1,xisz
        do j=1,zjsz
          do k=1,ny
!           Set the field per u(x,y,z) = x + z^2 + y^3
            u(k,j,i,n) = k*k*k + (zjst+j)*(zjst+j) + (xist+i)
          enddo
        enddo
      enddo
    enddo

  end subroutine initialize_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Decompose the global grid using the global problem size and the
!! number of MPI ranks.  Decomposition is always 2D on a Cartesian
!! topology.
!!
!! Subroutine sets xisz, zjsz, etc. with the decomposition results.
!!
!! \param[in] procnumber is the local MPI rank
!! \param[in] total      is the MPI communicator size
  subroutine initialize_problem(procnumber, total)

!   Input
    integer(4), intent(in) :: procnumber,total

!   Analogous to nx, ny, nz, and nx/2: artifacts of dealiasing
    integer(4), parameter :: nypad = ny
    integer(4), parameter :: nzpad = nz
    integer(4), parameter :: nxpad = nx, nxh = nx / 2

!   MPI process info
    integer, parameter :: MPI_PROC_NULL = -1
    integer, parameter :: MPI_STATUS_SIZE = 4
    integer ierr, dims(3),  cartid(3)
    logical periodic(3),remain_dims(3)
    integer iproc, jproc, ipid, jpid
    integer impid, ippid, jmpid, jppid

    integer, dimension(:), allocatable :: iist,iien,iisz
    integer, dimension(:), allocatable :: jjst,jjen,jjsz
    integer, dimension(:), allocatable :: kist,kien,kisz
    integer, dimension(:), allocatable :: kjst,kjen,kjsz

    integer mpi_comm_cart
    integer mpi_comm_row, mpi_comm_col

    logical iex

    dims(1) = 0
    dims(2) = 0

!   total is divided into a iproc x jproc stencil
    call MPI_Dims_create(total,2,dims,ierr)

    if(dims(1) .gt. dims(2)) then
      dims(1) = dims(2)
      dims(2) = total / dims(1)
    endif

    inquire(file='dims',exist=iex)
    if (iex) then
      if (procnumber.eq.0) write(6,*), 'Reading grid from file dims'
      open (999,file='dims')
      read (999,*) dims(1), dims(2)
      close (999)
    else
      if (procnumber.eq.0) write(6,*), 'Creating grid with mpi_dims_create'
    endif

    iproc = dims(1)
    jproc = dims(2)

    if (iproc*jproc.ne.total) then
      if (procnumber.eq.0) then
        write (6,*)  'ABORT: invalid choice of iproc x jproc!'
        write (6,*)  'Correct choices in the dims file'
        write (6,*) 'iproc,jproc,total=',iproc,jproc,total
        call MPI_ABORT(mpi_comm_world,1,ierr)
      end if
    end if

    periodic(1) = .false.
    periodic(2) = .false.

!   Create Cartesian processor grid and obtain new processor IDs
    call MPI_Cart_create(mpi_comm_world,2,dims,periodic, &
                         .false.,mpi_comm_cart,ierr)
    call MPI_Cart_coords(mpi_comm_cart,procnumber,2,cartid,ierr)
    ipid = cartid(1)
    jpid = cartid(2)

!   Compute neighbors in topology: i is east-west j is north-south
    impid = ipid - 1                            ! west neighbor
    ippid = ipid + 1                            ! east neighbor
    jmpid = jpid - 1                            ! north neighbor
    jppid = jpid + 1                            ! south neighbot
    if (ipid.eq.0)       impid = MPI_PROC_NULL  ! Boundaries
    if (jpid.eq.0)       jmpid = MPI_PROC_NULL
    if (ipid.eq.iproc-1) ippid = MPI_PROC_NULL
    if (jpid.eq.jproc-1) jppid = MPI_PROC_NULL

!   Find east-west row and north-south column communicators
    remain_dims(1) = .true.
    remain_dims(2) = .false.
    call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_row,ierr)
    remain_dims(1) = .false.
    remain_dims(2) = .true.
    call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_col,ierr)

!   Mapping i onto iproc, i onto jproc, j onto iproc etc.
    allocate (iist(0:iproc-1))
    allocate (iisz(0:iproc-1))
    allocate (iien(0:iproc-1))
    allocate (jjst(0:jproc-1))
    allocate (jjsz(0:jproc-1))
    allocate (jjen(0:jproc-1))
    allocate (kist(0:iproc-1))
    allocate (kisz(0:iproc-1))
    allocate (kien(0:iproc-1))
    allocate (kjst(0:jproc-1))
    allocate (kjsz(0:jproc-1))
    allocate (kjen(0:jproc-1))

!   Map 3-D data arrays onto 2-D process grid
!   (nx+2,ny+2,nz) => (iproc,jproc)
    call MapDataToProc(nxh,iproc,iist,iien,iisz)
    call MapDataToProc(nypad,jproc,jjst,jjen,jjsz)
    call MapDataToProc(nzpad,iproc,kist,kien,kisz)
    call MapDataToProc(nz,jproc,kjst,kjen,kjsz)

    xist = iist(ipid)   !x start
    xisz = iisz(ipid)   !x size
    xien = iien(ipid)   !x end

    yjst = jjst(jpid)
    yjsz = jjsz(jpid)
    yjen = jjen(jpid)

    zist = kist(ipid)
    zisz = kisz(ipid)
    zien = kien(ipid)

    zjsz = kjsz(jpid)
    zjen = kjen(jpid)
    zjst = kjst(jpid)

!   Finally, allocate local test field storage
    allocate (u(ny,zjsz,xisz,nc))

  end subroutine initialize_problem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine MapDataToProc (data,proc,st,en,sz)
    integer data,proc,st(0:proc-1),en(0:proc-1),sz(0:proc-1)
    integer i,size,nadd,size2
    size=data/proc
    nadd=mod(data,proc)
    size2=size
    if(nadd.ne.0) size2= size2+1
    st(0) = 1
    sz(0) = size2
    en(0) = size2
    if (proc .gt. 1) then
      do i=1,proc-1
        size2=size
        if (i.lt.nadd) size2=size2+1
        st(i) = st(i-1) + sz(i-1)
        sz(i) = size2
        en(i) = en(i-1) + size2
      enddo
      en(proc-1)= data
      sz(proc-1)= data-st(i-1)+1
    endif
  end subroutine MapDataToProc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module p3dfft_like_header
