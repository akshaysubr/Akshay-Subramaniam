program postprocess

    use globals, only: rkind
    use globals, only: parallel,lowmem,x_c,y_c,z_c,u,v,t1,tf
    use globals, only: nx,ny,nz,ax,ay,az,ax1,axn,ay1,ayn,az1,azn,xInd
    use mpi, only: comm,proc,nprocs,master
    use mpi, only: x1proc,xnproc,y1proc,ynproc,z1proc,znproc
    use interfaces, only: ddx,ddy
    implicit none
    include 'mpif.h'
    
    integer ierr
    integer irc
    integer status(MPI_STATUS_SIZE)

    real(kind=rkind) uder

    integer :: step

    ! Set the low memory flag to true
    parallel=.TRUE.

    !--Initialize MPI
    call MPI_INIT( ierr )
    comm = MPI_COMM_WORLD

    call MPI_COMM_RANK( comm, proc, ierr )
    call MPI_COMM_SIZE( comm, nprocs, ierr )

    ! Read input file, metadata file and setup globals
    call setup_postprocess

    if (  master ) then
        print*,"nx, ny, nz = ",nx,ny,nz
    end if
    
    ! print*,proc,": ax, ay, az = ",ax,ay,az
    ! print*,proc,": ax1, ay1, az1 = ",ax1,ay1,az1
    ! print*,proc,": axn, ayn, azn = ",axn,ayn,azn
    ! print*,proc,": x1proc, xnproc = ",x1proc,xnproc

    step = 0
    if (master) print*, "> Processing time step ", step

    call read_parallel_grid
    call read_parallel_data(step)

    ! if (proc == 1) print*,"Proc ",proc,": x_c(:,ay1,az1) = ",x_c(:,ay1,az1)
    ! if (master) print*,"Proc ",proc,": z_c(ax1,ay1,:) = ",z_c(ax1,ay1,:)
   
    ! Communicate data
    call CommunicateXBoundaryData
    call CommunicateYBoundaryData
    call CommunicateZBoundaryData

    ! Wait for communications to end
    call CommunicateXBoundaryWait
    call CommunicateYBoundaryWait
    call CommunicateZBoundaryWait

    ! if (proc == 1) print*,"Proc ",proc,": x_c(:,ay1,az1) = ",x_c(:,ay1,az1)
    ! if (master) print*,"Proc ",proc,": z_c(ax1,ay1,:) = ",z_c(ax1,ay1,:)

    u = 2.0D0*x_c + 1.0D0*y_c

    call ddx(u,v)
    if (master) print*,"Maximum % error in du/dx = ",1.0D2 * MAXVAL(ABS(2.0D0 - v(ax1:axn,ay1:ayn,az1:azn) ))/2.0D0

    call ddy(u,v)
    if (master) print*,"Maximum % error in du/dy = ",1.0D2 * MAXVAL(ABS(1.0D0 - v(ax1:axn,ay1:ayn,az1:azn) ))/1.0D0

    ! if (x1proc .and. y1proc .and. z1proc) then
    !     print*, "Wall velocities along z: "
    !     print*, v(axn,ay1,:)
    ! end if

    ! Loop through time steps
    ! do step=t1,tf
    !     print*, '> Processing time step ',step
    !     call mypostprocess(step)
    ! end do
    
    call cleanup_postprocess

    ! Finalize MPI
    call MPI_FINALIZE(irc)

end program postprocess
