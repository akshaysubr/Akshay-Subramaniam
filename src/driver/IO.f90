! Subroutine to read the input file and meta data file
! Also sets the global variables and allocates the iodata array
subroutine setup_postprocess

    use globals, only: px,py,pz,ppx,ppy,ppz,nprocs,procmap,invprocmap,lowmem
    use globals, only: nx,ny,nz,ns,nvars,ndim,xInd,yInd,zInd
    use globals, only: ax,ay,az,dx,dy,dz 
    use globals, only: t1,tf,dt,nsteps
    use globals, only: jobdir,flen,verbose
    use globals, only: iodata
    implicit none
    integer :: p,xp,yp,zp
    integer :: ioUnit
    integer :: i

    character(len=90) :: inputFile,procmapfile,plotmir
    character(len=90) :: dumchar

    NAMELIST /INPUT/ jobdir, ppx, ppy, ppz, t1, tf, dt, verbose

    ioUnit = 17

    IF (iargc() .EQ. 0) STOP 'Usage: ./readtest <input file>'
    CALL GETARG(1,inputFile)

    ! Read input file
    OPEN(UNIT=ioUnit, FILE=TRIM(inputFile), FORM='FORMATTED')
    READ(UNIT=ioUnit, NML=INPUT)
    if (verbose) print*,'JOB DIR: ',jobdir
    CLOSE(ioUnit)

    ! Read the processor map
    WRITE(procmapfile,*) TRIM(jobdir),'/procmap'
    OPEN(UNIT=ioUnit, FILE=TRIM(procmapfile), FORM='FORMATTED')
    READ(ioUnit,*) pz, py, px

    ! Skip header line
    READ(ioUnit,*)

    nprocs = px*py*pz

    allocate( procmap(nprocs,4) )
    allocate( invprocmap(px,py,pz) )
    DO i=1,nprocs
        READ(ioUnit,*) p,zp,yp,xp
        procmap(p+1,1) = p
        procmap(p+1,2) = xp
        procmap(p+1,3) = yp
        procmap(p+1,4) = zp
        invprocmap(xp+1,yp+1,zp+1) = p
    END DO
    CLOSE(ioUnit)

    ! Read plot.mir file to get # of variables and materials
    WRITE(plotmir,*) TRIM(jobdir),'/plot.mir'
    OPEN(UNIT=ioUnit, FILE=TRIM(plotmir), FORM='FORMATTED')

    ! Read metadata
    DO i=1,5; READ(ioUnit,*); END DO         ! Skip first 5 lines
    READ(ioUnit,*) dumchar,nx,ny,nz          ! Domain size
    DO i=1,2; READ(ioUnit,*); END DO         ! Skip next 2 lines
    READ(ioUnit,*) dumchar,dx,dy,dz          ! Grid spacing
    READ(ioUnit,*) dumchar, nvars            ! # of variables
    DO i=1,nvars; READ(ioUnit,*); END DO     ! Skip variable lines
    READ(ioUnit,*) dumchar, ns               ! # of species
    DO i=1,ns; READ(ioUnit,*); END DO        ! Skip material lines
    READ(ioUnit,*) dumchar, nsteps           ! # of timesteps

    nvars = nvars + 2 ! Since one of the vars is velocity

    CLOSE(ioUnit)

    call setup_globals
    if (lowmem) then
        call allocate_singleProc
    else
        call allocate_allProcs
    end if
    call setup_pointers

end subroutine setup_postprocess


! Subroutine to deallocate all the allocated global arrays
subroutine cleanup_postprocess

    use globals, only: iodata,procmap,invprocmap
    implicit none

    deallocate( iodata )
    deallocate( procmap )
    deallocate( invprocmap )

end subroutine cleanup_postprocess


! Read the data from a given processor and visualization directory
subroutine readProcData(vizdir,xp,yp,zp,procdata)

    use globals, only: nx,ny,nz,px,py,pz,ax,ay,az,nvars,ndim,invprocmap,flen,verbose
    use globals, only: xInd,yInd,zInd
    use globals, only: iodata
    implicit none
    integer, intent(in):: xp,yp,zp
    integer :: proc,pUnit,var
    character(len=flen), intent(in) :: vizdir
    character(len=flen) :: procfile
    real(kind=4), dimension(ax,ay,az,ndim-3), intent(out) :: procdata

    pUnit = 27

    ! Get processor ID
    proc = invprocmap(xp+1,yp+1,zp+1)

    ! Read in variables
    WRITE(procfile,'(2A,I6.6)') TRIM(vizdir),'/p',proc
    if(verbose) print*,'Reading file: ',TRIM(procfile)
    OPEN(UNIT=pUnit,FILE=TRIM(procfile),FORM='UNFORMATTED',STATUS='OLD')
    DO var=1,ndim-3
        READ(pUnit) procdata(:,:,:,var)
    END DO
    CLOSE(pUnit)

end subroutine readProcData


! Read the grid from a given processor
subroutine readProcGrid(xp,yp,zp,procgrid)

    use globals, only: nx,ny,nz,px,py,pz,ax,ay,az,nvars,ndim,invprocmap,jobdir,flen,verbose
    use globals, only: xInd,yInd,zInd
    use globals, only: iodata
    implicit none
    integer, intent(in):: xp,yp,zp
    integer :: proc,pUnit,var
    character(len=flen) :: procfile
    real(kind=4), dimension(ax,ay,az,3), intent(out) :: procgrid

    pUnit = 27

    ! Get processor ID
    proc = invprocmap(xp+1,yp+1,zp+1)

    ! Read in the grid
    WRITE(procfile,'(2A,I6.6)') TRIM(jobdir),'/grid/p',proc
    if (verbose) print*,'Reading file: ',TRIM(procfile)
    OPEN(UNIT=pUnit,FILE=TRIM(procfile),FORM='UNFORMATTED',STATUS='OLD')
    DO var=1,3
        READ(pUnit) procgrid(:,:,:,var)
    END DO
    CLOSE(pUnit)

end subroutine readProcGrid


! Read data from all procs and store in global iodata array
subroutine read_allProcs(step)
    use globals, only: nx,ny,nz,ax,ay,az,px,py,pz,xInd,yInd,zInd,ns,ndim,iodata
    use globals, only: jobdir,flen
    implicit none
    integer, intent(in) :: step
    integer :: xp,yp,zp
    character(len=flen) :: vizdir
    real(kind=4), dimension(:,:,:,:), allocatable :: procdata

    if (.not. allocated(iodata)) call allocate_allProcs
    allocate( procdata(ax,ay,az,ndim-3) )

    WRITE(vizdir,'(2A,I4.4)') TRIM(jobdir),'/vis',step

    do xp=0,px-1
        do yp=0,py-1
            do zp=0,pz-1
                call readProcData(vizdir,xp,yp,zp,procdata)
                iodata(xp*ax+1:(xp+1)*ax,yp*ay+1:(yp+1)*ay,zp*az+1:(zp+1)*az,1:ndim-3) = procdata
            end do
        end do
    end do

    deallocate( procdata )

end subroutine read_allProcs


! Read data from single proc and store in global iodata array
subroutine read_thisProc(step,xp,yp,zp)
    use globals, only: nx,ny,nz,ax,ay,az,px,py,pz,xInd,yInd,zInd,ns,ndim,iodata
    use globals, only: jobdir,flen
    implicit none
    integer, intent(in) :: step
    integer, intent(in) :: xp,yp,zp
    character(len=flen) :: vizdir
    real(kind=4), dimension(:,:,:,:), allocatable :: procdata

    if (.not. allocated(iodata)) call allocate_singleProc
    allocate( procdata(ax,ay,az,ndim-3) )

    WRITE(vizdir,'(2A,I4.4)') TRIM(jobdir),'/vis',step

    call readProcData(vizdir,xp,yp,zp,procdata)
    iodata(:,:,:,1:ndim-3) = procdata

    deallocate( procdata )

end subroutine read_thisProc


! Read grid from all the processors and store in global iodata array
subroutine read_allProcs_grid
    use globals, only: nx,ny,nz,ax,ay,az,px,py,pz,xInd,yInd,zInd,ns,ndim,iodata
    use globals, only: jobdir,flen
    implicit none
    integer :: xp,yp,zp
    real(kind=4), dimension(:,:,:,:), allocatable :: procgrid

    if (.not. allocated(iodata)) call allocate_allProcs
    allocate( procgrid(ax,ay,az,3) )

    do xp=0,px-1
        do yp=0,py-1
            do zp=0,pz-1
                call readProcGrid(xp,yp,zp,procgrid)
                iodata(xp*ax+1:(xp+1)*ax,yp*ay+1:(yp+1)*ay,zp*az+1:(zp+1)*az,ndim-2:ndim) = procgrid
            end do
        end do
    end do

    deallocate( procgrid )

end subroutine read_allProcs_grid


! Read grid from this processor and store in global iodata array
subroutine read_thisProc_grid(xp,yp,zp)
    use globals, only: nx,ny,nz,ax,ay,az,px,py,pz,xInd,yInd,zInd,ns,ndim,iodata
    use globals, only: jobdir,flen
    implicit none
    integer, intent(in) :: xp,yp,zp
    real(kind=4), dimension(:,:,:,:), allocatable :: procgrid

    if (.not. allocated(iodata)) call allocate_allProcs
    allocate( procgrid(ax,ay,az,3) )

    call readProcGrid(xp,yp,zp,procgrid)
    iodata(:,:,:,ndim-2:ndim) = procgrid

    deallocate( procgrid )

end subroutine read_thisProc_grid

subroutine read_parallel_data(step)
    use globals, only: nx,ny,nz,ax,ay,az,ax1,axn,ay1,ayn,az1,azn,nb,px,py,pz,ppx,ppy,ppz,ns,ndim,iodata,verbose,jobdir,flen,invprocmap
    use mpi, only: proc,xproc,yproc,zproc
    implicit none
    integer, intent(in) :: step
    
    integer :: p,xp,yp,zp
    integer :: x1,xn,y1,yn,z1,zn
    integer :: pUnit,var
    character(len=flen) :: vizdir
    character(len=flen) :: procfile

    pUnit = 27
    
    WRITE(vizdir,'(2A,I4.4)') TRIM(jobdir),'/vis',step

    do xp=xproc*ppx,(xproc+1)*ppx-1
        do yp=yproc*ppy,(yproc+1)*ppy-1
            do zp=zproc*ppz,(zproc+1)*ppz-1

                ! Get processor ID
                p = invprocmap(xp+1,yp+1,zp+1)

                x1 = ax1 + (xp-xproc*ppx)*nx/px
                xn = x1 + nx/px - 1
                y1 = ay1 + (yp-yproc*ppy)*ny/py
                yn = y1 + ny/py - 1
                z1 = az1 + (zp-zproc*ppz)*nz/pz
                zn = z1 + nz/pz - 1

                ! Read in variables
                WRITE(procfile,'(2A,I6.6)') TRIM(vizdir),'/p',p
                if(verbose) print*,'Reading file: ',TRIM(procfile)
                OPEN(UNIT=pUnit,FILE=TRIM(procfile),FORM='UNFORMATTED',STATUS='OLD')
                DO var=1,ndim-3
                    READ(pUnit) iodata(x1:xn,y1:yn,z1:zn,var)
                END DO
                CLOSE(pUnit)

            end do
        end do
    end do

end subroutine read_parallel_data

subroutine read_parallel_grid
    use globals, only: nx,ny,nz,ax,ay,az,ax1,axn,ay1,ayn,az1,azn,nb,px,py,pz,ppx,ppy,ppz,ns,ndim,iodata,verbose,jobdir,flen,invprocmap
    use mpi, only: xproc,yproc,zproc
    implicit none
    
    integer :: p,xp,yp,zp
    integer :: x1,xn,y1,yn,z1,zn
    integer :: pUnit,var
    character(len=flen) :: procfile

    pUnit = 27
    
    do xp=xproc*ppx,(xproc+1)*ppx-1
        do yp=yproc*ppy,(yproc+1)*ppy-1
            do zp=zproc*ppz,(zproc+1)*ppz-1

                ! Get processor ID
                p = invprocmap(xp+1,yp+1,zp+1)

                x1 = ax1 + (xp-xproc*ppx)*nx/px
                xn = x1 + nx/px - 1
                y1 = ay1 + (yp-yproc*ppy)*ny/py
                yn = y1 + ny/py - 1
                z1 = az1 + (zp-zproc*ppz)*nz/pz
                zn = z1 + nz/pz - 1

                ! Read in the grid
                WRITE(procfile,'(2A,I6.6)') TRIM(jobdir),'/grid/p',p
                if (verbose) print*,'Reading file: ',TRIM(procfile)
                OPEN(UNIT=pUnit,FILE=TRIM(procfile),FORM='UNFORMATTED',STATUS='OLD')
                DO var=1,3
                    READ(pUnit) iodata(x1:xn,y1:yn,z1:zn,ndim-3+var)
                END DO
                CLOSE(pUnit)

            end do
        end do
    end do

end subroutine read_parallel_grid
