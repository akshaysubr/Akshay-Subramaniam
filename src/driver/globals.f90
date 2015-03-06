module globals
    implicit none
    integer, parameter :: flen=120
    logical :: verbose=.FALSE.
    logical :: lowmem=.FALSE.
    logical :: parallel=.TRUE.

    integer :: px                                              ! No of procs in x
    integer :: py                                              ! No of procs in y
    integer :: pz                                              ! No of procs in z
    integer :: nprocs                                          ! Total no of procs

    integer :: ppx                                             ! No of viz x-procs per x proc 
    integer :: ppy                                             ! No of viz y-procs per y proc 
    integer :: ppz                                             ! No of viz z-procs per z proc 
    
    integer, parameter :: nb=3                                 ! No of boundary slices on each side (Now equal to maximum derivative required)
    integer :: nx                                              ! Total x points
    integer :: ny                                              ! Total y points
    integer :: nz                                              ! Total z points
    integer :: ax                                              ! Total x points per proc
    integer :: ay                                              ! Total y points per proc
    integer :: az                                              ! Total z points per proc
    integer :: ax1                                             ! First interior x point index per proc
    integer :: ay1                                             ! First interior y point index per proc
    integer :: az1                                             ! First interior z point index per proc
    integer :: axn                                             ! Last interior x point index per proc
    integer :: ayn                                             ! Last interior y point index per proc
    integer :: azn                                             ! Last interior z point index per proc
    double precision :: dx                                     ! Grid spacing in x direction
    double precision :: dy                                     ! Grid spacing in y direction
    double precision :: dz                                     ! Grid spacing in z direction

! ============== INPUT VARIABLES ====================
    character(len=flen) :: jobdir                              ! Job directory with all viz files 
    real(kind=4) :: dt                                         ! Time between viz dumps
    integer :: t1                                              ! First viz dump to be processed
    integer :: tf                                              ! Last viz dump to be processed
! ============== INPUT VARIABLES ====================
    
    integer :: nsteps                                          ! Total number of available viz dumps

    integer, dimension(:,:), allocatable :: procmap            ! Processor to grid mapping
    integer, dimension(:,:,:), allocatable :: invprocmap       ! Grid to processor mapping

    integer :: ndim                                            ! Total number of variables (including coordinates)
    integer :: nvars                                           ! Total number of flow variables
    integer :: ns                                              ! Number of species
    integer, parameter :: uInd=1                               ! x-velocity Index
    integer, parameter :: vInd=2                               ! y-velocity Index
    integer, parameter :: wInd=3                               ! z-velocity Index
    integer, parameter :: rhoInd=4                             ! Density Index
    integer, parameter :: eInd=5                               ! Energy Index
    integer, parameter :: pInd=6                               ! Pressure Index
    integer, parameter :: TInd=7                               ! Temp Index
    integer, parameter :: cInd=8                               ! Speed of sound Index
    integer, parameter :: muInd=9                              ! Shear visc Index
    integer, parameter :: bulkInd=10                           ! Bulk visc Index
    integer, parameter :: ktcInd=11                            ! Thermal cond Index
    integer, parameter :: DiffInd=12                           ! Species diffusion Index
    integer, parameter :: Y1Ind=13                             ! Species mass-fraction Index

    integer :: ncoords=3
    integer :: xInd                                            ! x-coordinate Index
    integer :: yInd                                            ! y-coordinate Index
    integer :: zInd                                            ! z-coordinate Index

    real(kind=4), dimension(:,:,:,:), allocatable, target :: iodata
    real(kind=4), dimension(:,:,:)  ,  pointer :: u             ! x velocity component
    real(kind=4), dimension(:,:,:)  ,  pointer :: v             ! y velocity component
    real(kind=4), dimension(:,:,:)  ,  pointer :: w             ! z velocity component
    real(kind=4), dimension(:,:,:)  ,  pointer :: rho           ! density
    real(kind=4), dimension(:,:,:)  ,  pointer :: e             ! internal energy
    real(kind=4), dimension(:,:,:)  ,  pointer :: p             ! pressure
    real(kind=4), dimension(:,:,:)  ,  pointer :: T             ! temperature
    real(kind=4), dimension(:,:,:)  ,  pointer :: c             ! speed of sound
    real(kind=4), dimension(:,:,:)  ,  pointer :: mu            ! shear viscosity
    real(kind=4), dimension(:,:,:)  ,  pointer :: bulk          ! bulk viscosity
    real(kind=4), dimension(:,:,:)  ,  pointer :: ktc           ! thermal conductivity
    real(kind=4), dimension(:,:,:)  ,  pointer :: Diff          ! diffusivity
    real(kind=4), dimension(:,:,:,:),  pointer :: Y             ! mass fractions
    real(kind=4), dimension(:,:,:)  ,  pointer :: x_c           ! x coordinate
    real(kind=4), dimension(:,:,:)  ,  pointer :: y_c           ! y coordinate
    real(kind=4), dimension(:,:,:)  ,  pointer :: z_c           ! z coordinate

    double precision, parameter :: pi=3.14159265358979323846d0,zero=0.0d0,one=1.0d0,two=2.0d0

end module globals



! Subroutine to setup all the global variables
subroutine setup_globals

    use globals, only: nx,ny,nz,nb,px,py,pz,ppx,ppy,ppz,t1,tf,nsteps
    use globals, only: ax,ay,az,ax1,axn,ay1,ayn,az1,azn
    use globals, only: xInd,yInd,zInd,ndim,nvars,ns,iodata
    use globals, only: parallel,lowmem
    use mpi, only: x1proc,xnproc,y1proc,ynproc,z1proc,znproc,bcount_x,bcount_y,bcount_z
    implicit none

    ! If parallel job, set the lowmem flag to true
    if (parallel) then
        lowmem = .TRUE.
        call setup_mpi
    end if

    ! Set points per processor
    ax = ppx * (nx/px) + 2*nb
    ay = ppy * (ny/py) + 2*nb
    az = ppz * (nz/pz) + 2*nb

    ax1 = nb+1
    axn = ax-nb
    ay1 = nb+1
    ayn = ay-nb
    az1 = nb+1
    azn = az-nb
   
    ! Case for boundary procs
    if (x1proc) then
        ax = ax - nb
        ax1 = 1
        axn = ax-nb
    end if
    if (xnproc) then
        ax = ax - nb
        axn = ax
    end if
    if (y1proc) then
        ay = ay - nb
        ay1 = 1
        ayn = ay-nb
    end if
    if (ynproc) then
        ay = ay - nb
        ayn = ay
    end if
    if (z1proc) then
        az = az - nb
        az1 = 1
        azn = az-nb
    end if
    if (znproc) then
        az = az - nb
        azn = az
    end if

    bcount_x = nb*ay*az
    bcount_y = nb*ax*az
    bcount_z = nb*ax*ay

    ! Set indices for coords
    xInd = nvars+ns+1
    yInd = xInd+1
    zInd = yInd+1

    ! Set total variable dimensionality
    ndim = zInd

    ! Check validity of processor arrangement. Only integer number of viz procs
    ! per compute proc allowed.
    if ( MOD(px,ppx) .NE. 0 ) then
        print*,'ERROR: Invalid number of viz x-procs',ppx,' per proc. Must be a multiple of total x procs, ',px
        stop
    else if ( MOD(py,ppy) .NE. 0 ) then
        print*,'ERROR: Invalid number of viz x-procs',ppx,' per proc. Must be a multiple of total x procs, ',px
        stop
    else if ( MOD(pz,ppz) .NE. 0 ) then
        print*,'ERROR: Invalid number of viz x-procs',ppx,' per proc. Must be a multiple of total x procs, ',px
        stop
    end if

    ! Check validity of start and end viz step values
    if (t1 .gt. tf) then
        print*,'ERROR: Invalid initial and final times (t1 > tf). Check input file and rerun'
        stop
    end if
    if (tf .gt. nsteps-1) then
        print*,'WARNING: Final time greater than max available time. Setting final time to last available time'
        tf = nsteps
    end if
    if (t1 .gt. nsteps-1) then
        print*,'ERROR: Invalid initial time t1 (t1 > max available time). Check input file and rerun'
        stop
    end if

end subroutine setup_globals



! Allocate memory for data from all procs
subroutine allocate_allProcs

    use globals, only: nx,ny,nz,ndim,iodata
    implicit none

    if (.not. allocated(iodata)) allocate( iodata(nx,ny,nz,ndim) )

end subroutine allocate_allProcs


! Allocate memory for data from single proc
subroutine allocate_singleProc

    use globals, only: nb,ax,ay,az,ndim,iodata,parallel
    use mpi, only: sendBuf_xl, sendBuf_yl, sendBuf_zl
    use mpi, only: sendBuf_xr, sendBuf_yr, sendBuf_zr
    use mpi, only: recvBuf_xl, recvBuf_yl, recvBuf_zl
    use mpi, only: recvBuf_xr, recvBuf_yr, recvBuf_zr
    implicit none

    if (.not. allocated(iodata)) allocate( iodata(ax,ay,az,ndim) )
    
    if (.not.allocated(sendBuf_xl)) allocate( sendBuf_xl(nb,ay,az,ndim) )
    if (.not.allocated(sendBuf_xr)) allocate( sendBuf_xr(nb,ay,az,ndim) )
    if (.not.allocated(sendBuf_yl)) allocate( sendBuf_yl(ax,nb,az,ndim) )
    if (.not.allocated(sendBuf_yr)) allocate( sendBuf_yr(ax,nb,az,ndim) )
    if (.not.allocated(sendBuf_zl)) allocate( sendBuf_zl(ax,ay,nb,ndim) )
    if (.not.allocated(sendBuf_zr)) allocate( sendBuf_zr(ax,ay,nb,ndim) )

    if (.not.allocated(recvBuf_xl)) allocate( recvBuf_xl(nb,ay,az,ndim) )
    if (.not.allocated(recvBuf_xr)) allocate( recvBuf_xr(nb,ay,az,ndim) )
    if (.not.allocated(recvBuf_yl)) allocate( recvBuf_yl(ax,nb,az,ndim) )
    if (.not.allocated(recvBuf_yr)) allocate( recvBuf_yr(ax,nb,az,ndim) )
    if (.not.allocated(recvBuf_zl)) allocate( recvBuf_zl(ax,ay,nb,ndim) )
    if (.not.allocated(recvBuf_zr)) allocate( recvBuf_zr(ax,ay,nb,ndim) )

end subroutine allocate_singleProc


! Setup pointers to iodata for easier data access
subroutine setup_pointers

    use globals, only: u,v,w,rho,e,p,T,c,mu,bulk,ktc,Diff,Y,x_c,y_c,z_c,iodata,ns
    use globals, only: uInd,vInd,wInd,rhoInd,eInd,pInd,TInd,cInd,muInd,bulkInd,ktcInd,DiffInd,Y1Ind,xInd,yInd,zInd
    implicit none

    u    => iodata(:,:,:,    uInd)
    v    => iodata(:,:,:,    vInd)
    w    => iodata(:,:,:,    wInd)
    rho  => iodata(:,:,:,  rhoInd)
    e    => iodata(:,:,:,    eInd)
    p    => iodata(:,:,:,    pInd)
    T    => iodata(:,:,:,    TInd)
    c    => iodata(:,:,:,    cInd)
    mu   => iodata(:,:,:,   muInd)
    bulk => iodata(:,:,:, bulkInd)
    ktc  => iodata(:,:,:,  ktcInd)
    Diff => iodata(:,:,:, DiffInd)
    Y    => iodata(:,:,:, Y1Ind:Y1Ind+ns)
    x_c  => iodata(:,:,:,    xInd)
    y_c  => iodata(:,:,:,    yInd)
    z_c  => iodata(:,:,:,    zInd)

end subroutine setup_pointers

subroutine setup_mpi
    use globals, only: nb,ax,ay,az,px,py,pz,ppx,ppy,ppz,procmap
    use mpi, only: proc,nprocs,xproc,yproc,zproc,npx,npy,npz
    use mpi, only: master,master_proc,x1proc,xnproc,y1proc,ynproc,z1proc,znproc
    use mpi, only: proc_xl,proc_xr,proc_yl,proc_yr,proc_zl,proc_zr,GetProcID
    implicit none
    integer :: xp,yp,zp

    npx = px/ppx
    npy = py/ppy
    npz = pz/ppz

    if ( npx*npy*npz .ne. nprocs ) then
        print*,'ERROR: Total number of procs not equal to product of procs in each direction. Check ppx,ppy,ppz values'
        print*,'ppx,ppy,ppz = ',ppx,ppy,ppz
        print*,' px, py, pz = ',px,py,pz
        stop
    end if

    xproc = MOD( proc,npx )
    yproc = MOD( proc/npx,npy )
    zproc = (proc/npx)/npy

    if (  proc == master_proc ) master = .TRUE.

    if ( xproc == 0 ) x1proc = .TRUE.
    if ( yproc == 0 ) y1proc = .TRUE.
    if ( zproc == 0 ) z1proc = .TRUE.

    if ( xproc == npx-1 ) xnproc = .TRUE.
    if ( yproc == npy-1 ) ynproc = .TRUE.
    if ( zproc == npz-1 ) znproc = .TRUE.

    if (.not. x1proc) proc_xl = GetProcID(xproc-1,yproc,zproc)
    if (.not. xnproc) proc_xr = GetProcID(xproc+1,yproc,zproc)
    if (.not. y1proc) proc_yl = GetProcID(xproc,yproc-1,zproc)
    if (.not. ynproc) proc_yr = GetProcID(xproc,yproc+1,zproc)
    if (.not. z1proc) proc_zl = GetProcID(xproc,yproc,zproc-1)
    if (.not. znproc) proc_zr = GetProcID(xproc,yproc,zproc+1)

end subroutine setup_mpi

subroutine CommunicateXBoundaryData
    use globals, only:iodata,ax,ay,az,ax1,axn,ay1,ayn,az1,azn,nb,ndim
    use mpi, only: proc,proc_xl,proc_xr,proc_yl,proc_yr,proc_zl,proc_zr,x1proc,xnproc,y1proc,ynproc,z1proc,znproc
    use mpi, only: comm,datatype,stat_size
    use mpi, only: reqSend_xl,reqSend_xr,reqSend_yl,reqSend_yr,reqSend_zl,reqSend_zr
    use mpi, only: reqRecv_xl,reqRecv_xr,reqRecv_yl,reqRecv_yr,reqRecv_zl,reqRecv_zr
    use mpi, only: sendBuf_xl, sendBuf_yl, sendBuf_zl
    use mpi, only: sendBuf_xr, sendBuf_yr, sendBuf_zr
    use mpi, only: recvBuf_xl, recvBuf_yl, recvBuf_zl
    use mpi, only: recvBuf_xr, recvBuf_yr, recvBuf_zr

    implicit none

    integer :: count,tag,ierr,dest,source
    integer :: statSend(stat_size), statRecv(stat_size)

    count = nb*ay*az*ndim
    tag = 1
    if (.not. xnproc) then
        ! Communicate right x boundary data
        dest = proc_xr
        ! print*,proc,': Sending iodata(',axn-nb+1,':',axn,')'
        sendBuf_xr = iodata(axn-nb+1:axn,:,:,:)
        ! print*,proc,': Copied send data to buffer'
        call MPI_ISEND(sendBuf_xr,count,datatype,dest,tag,comm,reqSend_xr,ierr)
        ! print*,'Sent right boundary data from proc ',proc,' to proc ',dest
    end if
    if (.not. x1proc) then
        ! Recieve left boundary data
        source = proc_xl
        ! print*,proc,': Receiving'
        call MPI_IRECV(recvBuf_xl,count,datatype,source,tag,comm,reqRecv_xl,ierr)
        ! print*,'Waiting for left boundary data from proc ',source,' to proc ',proc
    end if

    count = nb*ay*az*ndim
    tag = tag+1
    if (.not. x1proc) then
        ! Communicate left x boundary data
        dest = proc_xl
        sendBuf_xl = iodata(ax1:ax1+nb-1,:,:,:)
        call MPI_ISEND(sendBuf_xl,count,datatype,dest,tag,comm,reqSend_xl,ierr)
        ! print*,'Sent left boundary data from proc ',proc,' to proc ',dest
    end if
    if (.not. xnproc) then
        ! Recieve right x boundary data
        source = proc_xr
        call MPI_IRECV(recvBuf_xr,count,datatype,source,tag,comm,reqRecv_xr,ierr)
        ! print*,'Waiting for right boundary data from proc ',source,' to proc ',proc
    end if
 
    ! print*,proc,': At the end of CommunicateXBoundaryData'

end subroutine CommunicateXBoundaryData

subroutine CommunicateXBoundaryWait

    use globals, only: iodata,nb,ax1,axn,ax,ndim
    use mpi, only: stat_size
    use mpi, only: proc,x1proc,xnproc,y1proc,ynproc,z1proc,znproc
    use mpi, only: reqSend_xl,reqSend_xr,reqSend_yl,reqSend_yr,reqSend_zl,reqSend_zr
    use mpi, only: reqRecv_xl,reqRecv_xr,reqRecv_yl,reqRecv_yr,reqRecv_zl,reqRecv_zr
    use mpi, only: sendBuf_xl, sendBuf_yl, sendBuf_zl
    use mpi, only: sendBuf_xr, sendBuf_yr, sendBuf_zr
    use mpi, only: recvBuf_xl, recvBuf_yl, recvBuf_zl
    use mpi, only: recvBuf_xr, recvBuf_yr, recvBuf_zr
    
    implicit none

    integer :: ierr
    integer :: statSend_xl(stat_size), statRecv_xl(stat_size)
    integer :: statSend_xr(stat_size), statRecv_xr(stat_size)
    
    ! print*,'In CommunicateXBoundaryWait'

    if (.not. xnproc) call MPI_WAIT(reqSend_xr,statSend_xr,ierr)
    if (.not. x1proc) then
        call MPI_WAIT(reqRecv_xl,statRecv_xl,ierr)
        ! print*,proc,': Recieved buffer. Now copying data to iodata'
        iodata(1:nb,:,:,:) = recvBuf_xl
    end if

    if (.not. x1proc) call MPI_WAIT(reqSend_xl,statSend_xl,ierr)
    if (.not. xnproc) then
        call MPI_WAIT(reqRecv_xr,statRecv_xr,ierr)
        ! print*,proc,': Recieved buffer. Now copying data to iodata'
        iodata(axn+1:ax,:,:,:) = recvBuf_xr
    end if

end subroutine CommunicateXBoundaryWait

subroutine CommunicateYBoundaryData
    use globals, only:iodata,ax,ay,az,ax1,axn,ay1,ayn,az1,azn,nb,ndim
    use mpi, only: proc,proc_xl,proc_xr,proc_yl,proc_yr,proc_zl,proc_zr,x1proc,xnproc,y1proc,ynproc,z1proc,znproc
    use mpi, only: comm,datatype,stat_size
    use mpi, only: reqSend_xl,reqSend_xr,reqSend_yl,reqSend_yr,reqSend_zl,reqSend_zr
    use mpi, only: reqRecv_xl,reqRecv_xr,reqRecv_yl,reqRecv_yr,reqRecv_zl,reqRecv_zr
    use mpi, only: sendBuf_xl, sendBuf_yl, sendBuf_zl
    use mpi, only: sendBuf_xr, sendBuf_yr, sendBuf_zr
    use mpi, only: recvBuf_xl, recvBuf_yl, recvBuf_zl
    use mpi, only: recvBuf_xr, recvBuf_yr, recvBuf_zr

    implicit none

    integer :: count,tag,ierr,dest,source
    integer :: statSend(stat_size), statRecv(stat_size)

    count = nb*ax*az*ndim
    tag = 11
    if (.not. ynproc) then
        ! Communicate right x boundary data
        dest = proc_yr
        ! print*,proc,': Sending iodata(',ayn-nb+1,':',ayn,')'
        sendBuf_yr = iodata(:,ayn-nb+1:ayn,:,:)
        ! print*,proc,': Copied send data to buffer'
        call MPI_ISEND(sendBuf_yr,count,datatype,dest,tag,comm,reqSend_yr,ierr)
        ! print*,'Sent right boundary data from proc ',proc,' to proc ',dest
    end if
    if (.not. y1proc) then
        ! Recieve left boundary data
        source = proc_yl
        ! print*,proc,': Receiving'
        call MPI_IRECV(recvBuf_yl,count,datatype,source,tag,comm,reqRecv_yl,ierr)
        ! print*,'Waiting for left boundary data from proc ',source,' to proc ',proc
    end if

    count = nb*ax*az*ndim
    tag = tag+1
    if (.not. y1proc) then
        ! Communicate left x boundary data
        dest = proc_yl
        sendBuf_yl = iodata(:,ay1:ay1+nb-1,:,:)
        call MPI_ISEND(sendBuf_yl,count,datatype,dest,tag,comm,reqSend_yl,ierr)
        ! print*,'Sent left boundary data from proc ',proc,' to proc ',dest
    end if
    if (.not. ynproc) then
        ! Recieve right x boundary data
        source = proc_yr
        call MPI_IRECV(recvBuf_yr,count,datatype,source,tag,comm,reqRecv_yr,ierr)
        ! print*,'Waiting for right boundary data from proc ',source,' to proc ',proc
    end if

    ! print*,proc,': At the end of CommunicateYBoundaryData'

end subroutine CommunicateYBoundaryData

subroutine CommunicateYBoundaryWait

    use globals, only: iodata,nb,ay1,ayn,ay,ndim
    use mpi, only: stat_size
    use mpi, only: proc,x1proc,xnproc,y1proc,ynproc,z1proc,znproc
    use mpi, only: reqSend_xl,reqSend_xr,reqSend_yl,reqSend_yr,reqSend_zl,reqSend_zr
    use mpi, only: reqRecv_xl,reqRecv_xr,reqRecv_yl,reqRecv_yr,reqRecv_zl,reqRecv_zr
    use mpi, only: sendBuf_xl, sendBuf_yl, sendBuf_zl
    use mpi, only: sendBuf_xr, sendBuf_yr, sendBuf_zr
    use mpi, only: recvBuf_xl, recvBuf_yl, recvBuf_zl
    use mpi, only: recvBuf_xr, recvBuf_yr, recvBuf_zr
    
    implicit none

    integer :: ierr
    integer :: statSend_yl(stat_size), statRecv_yl(stat_size)
    integer :: statSend_yr(stat_size), statRecv_yr(stat_size)
    
    ! print*,'In CommunicateYBoundaryWait'

    if (.not. ynproc) call MPI_WAIT(reqSend_yr,statSend_yr,ierr)
    if (.not. y1proc) then
        call MPI_WAIT(reqRecv_yl,statRecv_yl,ierr)
        ! print*,proc,': Recieved buffer. Now copying data to iodata'
        iodata(:,1:nb,:,:) = recvBuf_yl
    end if

    if (.not. y1proc) call MPI_WAIT(reqSend_yl,statSend_yl,ierr)
    if (.not. ynproc) then
        call MPI_WAIT(reqRecv_yr,statRecv_yr,ierr)
        ! print*,proc,': Recieved buffer. Now copying data to iodata'
        iodata(:,ayn+1:ay,:,:) = recvBuf_yr
    end if

end subroutine CommunicateYBoundaryWait

subroutine CommunicateZBoundaryData
    use globals, only:iodata,ax,ay,az,ax1,axn,ay1,ayn,az1,azn,nb,ndim
    use mpi, only: proc,proc_xl,proc_xr,proc_yl,proc_yr,proc_zl,proc_zr,x1proc,xnproc,y1proc,ynproc,z1proc,znproc
    use mpi, only: comm,datatype,stat_size
    use mpi, only: reqSend_xl,reqSend_xr,reqSend_yl,reqSend_yr,reqSend_zl,reqSend_zr
    use mpi, only: reqRecv_xl,reqRecv_xr,reqRecv_yl,reqRecv_yr,reqRecv_zl,reqRecv_zr
    use mpi, only: sendBuf_xl, sendBuf_yl, sendBuf_zl
    use mpi, only: sendBuf_xr, sendBuf_yr, sendBuf_zr
    use mpi, only: recvBuf_xl, recvBuf_yl, recvBuf_zl
    use mpi, only: recvBuf_xr, recvBuf_yr, recvBuf_zr

    implicit none

    integer :: count,tag,ierr,dest,source
    integer :: statSend(stat_size), statRecv(stat_size)

    count = nb*ax*ay*ndim
    tag = 21
    if (.not. znproc) then
        ! Communicate right x boundary data
        dest = proc_zr
        ! print*,proc,': Sending iodata(',azn-nb+1,':',azn,')'
        sendBuf_zr = iodata(:,:,azn-nb+1:azn,:)
        ! print*,proc,': Copied send data to buffer'
        call MPI_ISEND(sendBuf_zr,count,datatype,dest,tag,comm,reqSend_zr,ierr)
        ! print*,'Sent right boundary data from proc ',proc,' to proc ',dest
    end if
    if (.not. z1proc) then
        ! Recieve left boundary data
        source = proc_zl
        ! print*,proc,': Receiving'
        call MPI_IRECV(recvBuf_zl,count,datatype,source,tag,comm,reqRecv_zl,ierr)
        ! print*,'Waiting for left boundary data from proc ',source,' to proc ',proc
    end if

    count = nb*ax*ay*ndim
    tag = tag+1
    if (.not. z1proc) then
        ! Communicate left x boundary data
        dest = proc_zl
        sendBuf_zl = iodata(:,:,az1:az1+nb-1,:)
        call MPI_ISEND(sendBuf_zl,count,datatype,dest,tag,comm,reqSend_zl,ierr)
        ! print*,'Sent left boundary data from proc ',proc,' to proc ',dest
    end if
    if (.not. znproc) then
        ! Recieve right x boundary data
        source = proc_zr
        call MPI_IRECV(recvBuf_zr,count,datatype,source,tag,comm,reqRecv_zr,ierr)
        ! print*,'Waiting for right boundary data from proc ',source,' to proc ',proc
    end if

    ! print*,proc,': At the end of CommunicateZBoundaryData'

end subroutine CommunicateZBoundaryData

subroutine CommunicateZBoundaryWait

    use globals, only: iodata,nb,az1,azn,az,ndim
    use mpi, only: stat_size
    use mpi, only: proc,x1proc,xnproc,y1proc,ynproc,z1proc,znproc
    use mpi, only: reqSend_xl,reqSend_xr,reqSend_yl,reqSend_yr,reqSend_zl,reqSend_zr
    use mpi, only: reqRecv_xl,reqRecv_xr,reqRecv_yl,reqRecv_yr,reqRecv_zl,reqRecv_zr
    use mpi, only: sendBuf_xl, sendBuf_yl, sendBuf_zl
    use mpi, only: sendBuf_xr, sendBuf_yr, sendBuf_zr
    use mpi, only: recvBuf_xl, recvBuf_yl, recvBuf_zl
    use mpi, only: recvBuf_xr, recvBuf_yr, recvBuf_zr
    
    implicit none

    integer :: ierr
    integer :: statSend_zl(stat_size), statRecv_zl(stat_size)
    integer :: statSend_zr(stat_size), statRecv_zr(stat_size)
    
    ! print*,'In CommunicateZBoundaryWait'

    if (.not. znproc) call MPI_WAIT(reqSend_zr,statSend_zr,ierr)
    if (.not. z1proc) then
        call MPI_WAIT(reqRecv_zl,statRecv_zl,ierr)
        ! print*,proc,': Recieved buffer. Now copying data to iodata'
        iodata(:,:,1:nb,:) = recvBuf_zl
    end if

    if (.not. z1proc) call MPI_WAIT(reqSend_zl,statSend_zl,ierr)
    if (.not. znproc) then
        call MPI_WAIT(reqRecv_zr,statRecv_zr,ierr)
        ! print*,proc,': Recieved buffer. Now copying data to iodata'
        iodata(:,:,azn+1:az,:) = recvBuf_zr
    end if

end subroutine CommunicateZBoundaryWait
