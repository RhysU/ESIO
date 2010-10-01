!#
!# Module that contains the problem in use
!#
!# problem is a box with lengths ny,nx,nz
!# the indicies actually go from 1 -> ny+1, etc.
!#
!# I doubt that one needs to muck around much here
!# The variables you need: u (field), xist,xien,zjst,zjen (coordinates of field)
!#
module header
  use mpi
  implicit none ! this setting is global to module: applies to each subroutine

  ! 5 megs write size
  integer(4), parameter :: nx=50     ! x,y,z lengths (these wont change)
  integer(4), parameter :: ny=1000   ! these lengths are for physical space!
  integer(4), parameter :: nz=25

  integer(4), parameter :: nc=1      ! number of velocity components

  integer(4) :: xist,xjst,xisz,xjsz,xien,xjen ! pencil widths and heights
  integer(4) :: yist,yjst,yisz,yjsz,yien,yjen
  integer(4) :: zist,zjst,zisz,zjsz,zien,zjen

  real(8),dimension(:,:,:,:),allocatable :: u ! test velocity field: u(ny,zjsz,xisz,nc)
  real(8),dimension(:,:,:,:),allocatable :: bigu ! full velocity field: u(ny,nz,nx,nc)

  !public (anyone using the module gets these)
  public  :: u, xisz, zjsz, ny, nc                ! velocity field data - local
  public  :: xist,xien,zjst,zjen                  ! absolute coordinates of field
  public  :: nz,nx                                ! size of field
  public  :: initialize_problem, initialize_field ! public subroutines

contains

  !##################################################
  !#
  !#    Subroutine initialize_field
  !#
  !# This routine should generate a local
  !# initial field for each processor
  !#
  !# For reference, xisz (x-size) -- local coordinates
  !#                zjsz (z-size) -- on each processor
  !#                ny
  !#
  !# ny,nx,nz are the total sizes
  !#
  !#           and,
  !#                xist xien    (x start -- x finish)
  !#                zjst zjen    -- absolute coordinates
  !#                ny           -- on each proc
  !#
  !#  Lets set the field as follows:
  !#  f(x,y,z) = x + z^2 + y^3
  !#
  !##################################################
  subroutine initialize_field(proc)

    integer(4) :: i,j,k,n,proc

    !this loop is stride-1: u(ny,zjsz,xisz)
    !write(6,*) 'xisz, zjsz, ny: ', xisz, zjsz, ny

    do n=1,nc
       do i=1,xisz
          do j=1,zjsz
             do k=1,ny
!                u(k,j,i,n) = 1
                u(k,j,i,n) = k*k*k + (zjst+j)*(zjst+j) + (xist+i)
             enddo
          enddo
       enddo
    enddo

  end subroutine initialize_field

  !##################################################
  !#
  !# Subroutine initialize_problem
  !#
  !# this subroutine breaks the grid up based on
  !# how large the total problem is, as well as the
  !# number of processors being used
  !# we are always using a 2d-decomposition on
  !# a rectangular grid
  !#
  !# subroutine should set xisz, zjsz, etc.
  !#
  !# INPUT: procnumber is individual process number
  !#        total is number of processors
  !#
  !# Note: this routine has lots of junk in it
  !#
  !##################################################
  subroutine initialize_problem(procnumber,total)

      ! input
      integer(4), intent(in) :: procnumber,total

      !dont worry about these --they are just like ny,nz,nx, and nx/2
      integer(4) :: nxpad, nypad, nzpad,nxh

      ! mpi process info
      integer, parameter :: MPI_PROC_NULL = -1
      integer, parameter :: MPI_STATUS_SIZE = 4
!      integer :: mpi_comm_world
      integer ierr, dims(3),  cartid(3)
      logical periodic(3),remain_dims(3)
      integer mpierr
      integer iproc, jproc, ipid, jpid
      integer impid, ippid, jmpid, jppid

      real(8) mpixx

      integer, dimension(:), allocatable :: iist,iien,iisz
      integer, dimension(:), allocatable :: jjst,jjen,jjsz
      integer, dimension(:), allocatable :: kist,kien,kisz
      integer, dimension(:), allocatable :: kistpad,kienpad,kiszpad
      integer, dimension(:), allocatable :: kjst,kjen,kjsz

      !array of one's
      integer, dimension(:), allocatable :: oneAr

      integer :: padx,pady,padz,szmax,szf
      integer mpi_comm_cart
      integer mpi_comm_row, mpi_comm_col

      ! mpi derived data types for implementing alltoallv using send-recvs
      ! probably dont need any of this
      integer IfCntMax,KfCntMax,JrCntMax,KrCntMax
      logical IfCntUneven,KfCntUneven,JrCntUneven,KrCntUneven
      integer,dimension(:),allocatable:: IfSndCnts,IfSndStrt
      integer,dimension(:),allocatable:: IfRcvCnts,IfRcvStrt
      integer,dimension(:),allocatable:: KfSndCnts,KfSndStrt
      integer,dimension(:),allocatable:: KfRcvCnts,KfRcvStrt
      integer,dimension(:),allocatable:: JrSndCnts,JrSndStrt
      integer,dimension(:),allocatable:: JrRcvCnts,JrRcvStrt
      integer,dimension(:),allocatable:: KrSndCnts,KrSndStrt
      integer,dimension(:),allocatable:: KrRcvCnts,KrRcvStrt
      integer,dimension(:,:),allocatable:: status

      integer i,j,k,n1
      logical iex

      nypad=ny  ! dont worry about why we are using different vars for these
      nzpad=nz  ! (related to dealiasing)
      nxpad=nx
      nxh=nx/2  ! ; nxhp=nxh+1 ; nxhp2=2*nxhp

      allocate(status(MPI_STATUS_SIZE,total))

      dims(1) = 0
      dims(2) = 0

      !    total is divided into a iproc x jproc stencil
      !
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
!            call MPI_ABORT(mpi_comm_world,ierr)
            !some mysterious error here, possibly related to mpi_comm_world?
         end if
      end if

      periodic(1) = .false.
      periodic(2) = .false.

      ! creating cartesian processor grid
      !    MPI_CART_CREATE(mpi_comm_world,int ndims, int *dims,int *periods
      call MPI_Cart_create(mpi_comm_world,2,dims,periodic,.false.,mpi_comm_cart,ierr)

      ! Obtaining process ids with in the cartesian grid
      call MPI_Cart_coords(mpi_comm_cart,procnumber,2,cartid,ierr)

      ! process with a linear id of 5 may have cartid of (3,1)
      ipid = cartid(1)
      jpid = cartid(2)

      ! here i is east-west j is north-south
      ! impid is west neighbour ippid is east neighbour and so on
      impid = ipid - 1
      ippid = ipid + 1
      jmpid = jpid - 1
      jppid = jpid + 1
      !boundary processes
      if (ipid.eq.0) impid = MPI_PROC_NULL
      if (jpid.eq.0) jmpid = MPI_PROC_NULL
      if (ipid.eq.iproc-1) ippid = MPI_PROC_NULL
      if (jpid.eq.jproc-1) jppid = MPI_PROC_NULL
      ! using cart comworld create east-west(row) sub comworld

      remain_dims(1) = .true.
      remain_dims(2) = .false.
      call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_row,ierr)
      ! using cart comworld create north-south(column) sub comworld
      remain_dims(1) = .false.
      remain_dims(2) = .true.
      call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_col,ierr)

      ! mapping i onto iproc, i onto jproc, j onto iproc etc.
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

      !array of one's
      allocate (oneAr(0:total))
      do i=0,total
         oneAr(i) = 1
      enddo
      !
      !Mapping 3-D data arrays onto 2-D process grid
      ! (nx+2,ny+2,nz) => (iproc,jproc)
      !
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

      !write(6,*) "xist: ", xist
      !write(6,*) "xisz: ", xisz
      !write(6,*) "xien: ", xien

      !write(6,*) "yjst: ", yjst
      !write(6,*) "yjsz: ", yjsz
      !write(6,*) "yjen: ", yjen

      !write(6,*) "zist: ", zist
      !write(6,*) "zisz: ", zisz
      !write(6,*) "zien: ", zien

      !write(6,*) "zjst: ", zjst
      !write(6,*) "zjsz: ", zjsz
      !write(6,*) "zjen: ", zjen

      !      if(iproc .eq. 1) then
      !         xisz = nxhp
      !      endif

      !  mpi derived data types

      allocate (IfSndCnts(0:iproc-1))
      allocate (IfRcvCnts(0:iproc-1))
      allocate (KfSndCnts(0:jproc-1))
      allocate (KfRcvCnts(0:jproc-1))

      allocate (JrSndCnts(0:jproc-1))
      allocate (JrRcvCnts(0:jproc-1))

      allocate (KrSndCnts(0:iproc-1))
      allocate (KrRcvCnts(0:iproc-1))

      allocate (IfSndStrt(0:iproc-1))
      allocate (IfRcvStrt(0:iproc-1))
      allocate (KfSndStrt(0:jproc-1))
      allocate (KfRcvStrt(0:jproc-1))

      allocate (JrSndStrt(0:jproc-1))
      allocate (JrRcvStrt(0:jproc-1))

      allocate (KrSndStrt(0:iproc-1))
      allocate (KrRcvStrt(0:iproc-1))

      ! now, allocate velocity field
      allocate (u(ny,zjsz,xisz,nc))

  end subroutine initialize_problem


  !##################################################
  !# Subroutine MapDataToProc
  !#
  !# called by initialize_problem
  !#
  !#
  !##################################################

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
    !
  end subroutine MapDataToProc

  !##################################################
  !#    Subroutine buildbig
  !#
  !# called by main
  !# test program that builds the field pieces into
  !# one big field
  !#
  !#  Going to require some MPI calls
  !#
  !#  DONT USE for huge fields --- memory req too high
  !#
  !##################################################

  subroutine buildbig

    allocate (bigu(ny,nz,nx,nc))

  end subroutine buildbig
  !##################################################
  !#    Subroutine buildbig
  !#    same as build but not allocating
  !#    called by main
  !##################################################
  subroutine updatebig


  end subroutine updatebig

  end module header
