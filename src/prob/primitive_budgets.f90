program postprocess

    use globals, only: parallel,lowmem,x_c,y_c,z_c,t1,tf
    use globals, only: nx,ny,nz,ax,ay,az,ax1,axn,ay1,ayn,az1,azn,xInd
    use mpi, only: comm,proc,nprocs,master
    use mpi, only: x1proc,xnproc
    implicit none
    include 'mpif.h'
    
    integer ierr
    integer irc
    integer status(MPI_STATUS_SIZE)

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
        print*,"Grid size:"
        print*,"  nx, ny, nz = ",nx,ny,nz
        print*,"Total number of processors = ",nprocs
        print*,""
    end if
    
    ! Loop through time steps
    do step=t1,tf

        if (master) print*, '> Processing time step ',step
    
        call read_parallel_grid
        call read_parallel_data(step)

        ! Communicate data
        call CommunicateXBoundaryData
        call CommunicateYBoundaryData
        call CommunicateZBoundaryData

        !!! Need to add architecture to do some computations here !!!

        ! Wait for communications to end
        call CommunicateXBoundaryWait
        call CommunicateYBoundaryWait
        call CommunicateZBoundaryWait

        ! Now call the problem-specific postprocess routine
        call mypostprocess(step)
    end do
    
    call cleanup_postprocess

    ! Finalize MPI
    call MPI_FINALIZE(irc)

end program postprocess

subroutine mypostprocess(step)

    use globals, only: u,v,w,rho,p,e,T,mu,bulk,ktc,x_c,y_c,z_c,iodata
    use globals, only: dx,dy,dz,nx,ny,nz,px,py,pz,ax,ay,az
    use globals, only: t1,tf,dt,flen,jobdir
    use globals, only: ax1,ay1,az1,axn,ayn,azn
    use mpi, only: proc,master,master_proc,comm
    use interfaces, only: ddx,ddy,ddz,integrate,grad,curl,crossprod
    use functions, only: calc_vorticity
    
    implicit none
    
    include 'mpif.h'
    
    integer, intent(in) :: step
    integer :: xp,yp,zp,i,ierr
    real(kind=4), dimension(:,:,:), allocatable :: Em,dilatation,tmp
    real(kind=4), dimension(:,:,:,:), allocatable :: uder,vder,wder,Tder
    real(kind=4), dimension(:,:,:,:), allocatable :: x_flux,y_flux,z_flux,rhs
    real(kind=4), dimension(5) :: rhs_intg, rhs_intg_total
    character(len=flen) :: statsfile
    integer, parameter :: statsUnit=37
    logical, parameter :: writetofile=.FALSE.

    if (master) then
        if (writetofile) then
            WRITE(statsfile,'(2A)') TRIM(jobdir),'/primitive_budget.dat'
            if (step == t1) then
                OPEN(UNIT=statsUnit,FILE=statsfile,FORM='FORMATTED',STATUS='REPLACE')
            else
                OPEN(UNIT=statsUnit,FILE=statsfile,FORM='FORMATTED',POSITION='APPEND')
            end if
        end if
    end if

    ! Allocate required intermediate arrays
    allocate(     x_flux( SIZE(u,1), SIZE(u,2), SiZE(u,3), 5 ) )
    allocate(     y_flux( SIZE(u,1), SIZE(u,2), SiZE(u,3), 5 ) )
    allocate(     z_flux( SIZE(u,1), SIZE(u,2), SiZE(u,3), 5 ) )
    allocate(        rhs( SIZE(u,1), SIZE(u,2), SiZE(u,3), 5 ) )
    allocate(       uder( SIZE(u,1), SIZE(u,2), SiZE(u,3), 3 ) )
    allocate(       vder( SIZE(u,1), SIZE(u,2), SiZE(u,3), 3 ) )
    allocate(       wder( SIZE(u,1), SIZE(u,2), SiZE(u,3), 3 ) )
    allocate(       Tder( SIZE(u,1), SIZE(u,2), SiZE(u,3), 3 ) )
    allocate(         Em( SIZE(u,1), SIZE(u,2), SiZE(u,3)    ) )
    allocate( dilatation( SIZE(u,1), SIZE(u,2), SiZE(u,3)    ) )
    allocate(        tmp( SIZE(u,1), SIZE(u,2), SiZE(u,3)    ) )

    ! Calc fluxes

    Em = rho*(e + 0.5D0*u*u + 0.5D0*v*v + 0.5D0*w*w)

    x_flux(:,:,:,1) = rho*u
    x_flux(:,:,:,2) = rho*u*u + p
    x_flux(:,:,:,3) = rho*u*v
    x_flux(:,:,:,4) = rho*u*w
    x_flux(:,:,:,5) = (Em + p)*u

    y_flux(:,:,:,1) = rho*v
    y_flux(:,:,:,2) = rho*v*u
    y_flux(:,:,:,3) = rho*v*v + p
    y_flux(:,:,:,4) = rho*v*w
    y_flux(:,:,:,5) = (Em + p)*v

    z_flux(:,:,:,1) = rho*w
    z_flux(:,:,:,2) = rho*w*u
    z_flux(:,:,:,3) = rho*w*v
    z_flux(:,:,:,4) = rho*w*w + p
    z_flux(:,:,:,5) = (Em + p)*w

    ! Calc Sij and dilatation
    call ddx(u,uder(:,:,:,1))
    call ddy(u,uder(:,:,:,2))
    call ddz(u,uder(:,:,:,3))

    call ddx(v,vder(:,:,:,1))
    call ddy(v,vder(:,:,:,2))
    call ddz(v,vder(:,:,:,3))

    call ddx(w,wder(:,:,:,1))
    call ddy(w,wder(:,:,:,2))
    call ddz(w,wder(:,:,:,3))

    ! Get gradient of T
    call ddx(T,Tder(:,:,:,1))
    call ddy(T,Tder(:,:,:,2))
    call ddz(T,Tder(:,:,:,3))

    ! Get heat flux in Tder
    Tder(:,:,:,1) = -ktc*Tder(:,:,:,1)
    Tder(:,:,:,2) = -ktc*Tder(:,:,:,2)
    Tder(:,:,:,3) = -ktc*Tder(:,:,:,3)

    ! Get the Sij tensor in uder, vder and wder
    uder(:,:,:,2) = 0.5D0*(uder(:,:,:,2)+vder(:,:,:,1))
    vder(:,:,:,1) = uder(:,:,:,2)
    uder(:,:,:,3) = 0.5D0*(uder(:,:,:,3)+wder(:,:,:,1))
    wder(:,:,:,1) = uder(:,:,:,3)
    vder(:,:,:,3) = 0.5D0*(vder(:,:,:,3)+wder(:,:,:,2))
    wder(:,:,:,2) = vder(:,:,:,3)

    dilatation = uder(:,:,:,1)+vder(:,:,:,2)+wder(:,:,:,3)

    ! Now get the viscous stress tensor (Put rows in uder, vder, wder)
    do i=1,3
        uder(:,:,:,i) = 2.0D0*mu*uder(:,:,:,i) 
        vder(:,:,:,i) = 2.0D0*mu*vder(:,:,:,i) 
        wder(:,:,:,i) = 2.0D0*mu*wder(:,:,:,i) 
    end do
    uder(:,:,:,1) = uder(:,:,:,1) - (2.0D0*mu/3.0D0 - bulk)*dilatation
    vder(:,:,:,2) = vder(:,:,:,2) - (2.0D0*mu/3.0D0 - bulk)*dilatation
    wder(:,:,:,3) = wder(:,:,:,3) - (2.0D0*mu/3.0D0 - bulk)*dilatation

    ! Update fluxes with viscous terms

    x_flux(:,:,:,2) = x_flux(:,:,:,2) - uder(:,:,:,1)
    x_flux(:,:,:,3) = x_flux(:,:,:,3) - uder(:,:,:,2)
    x_flux(:,:,:,4) = x_flux(:,:,:,4) - uder(:,:,:,3)
    x_flux(:,:,:,5) = x_flux(:,:,:,5) - u*uder(:,:,:,1) - v*uder(:,:,:,2) - w*uder(:,:,:,3) - Tder(:,:,:,1)

    y_flux(:,:,:,2) = y_flux(:,:,:,2) - vder(:,:,:,1)
    y_flux(:,:,:,3) = y_flux(:,:,:,3) - vder(:,:,:,2)
    y_flux(:,:,:,4) = y_flux(:,:,:,4) - vder(:,:,:,3)
    y_flux(:,:,:,5) = y_flux(:,:,:,5) - u*vder(:,:,:,1) - v*vder(:,:,:,2) - w*vder(:,:,:,3) - Tder(:,:,:,2)

    z_flux(:,:,:,2) = x_flux(:,:,:,2) - uder(:,:,:,1)
    z_flux(:,:,:,3) = x_flux(:,:,:,3) - uder(:,:,:,2)
    z_flux(:,:,:,4) = x_flux(:,:,:,4) - uder(:,:,:,3)
    z_flux(:,:,:,5) = x_flux(:,:,:,5) - u*wder(:,:,:,1) - v*wder(:,:,:,2) - w*wder(:,:,:,3) - Tder(:,:,:,3)

    ! Get RHS = - div(fluxes)
    do i=1,5
        call ddx(x_flux(:,:,:,i),rhs(:,:,:,i))
        rhs(:,:,:,i) = -rhs(:,:,:,i)

        call ddy(y_flux(:,:,:,i),tmp)
        rhs(:,:,:,i) = rhs(:,:,:,i) - tmp

        call ddz(z_flux(:,:,:,i),tmp)
        rhs(:,:,:,i) = rhs(:,:,:,i) - tmp

        rhs_intg(i) = integrate(rhs(:,:,:,i))
    end do

    ! Now, reduce from all MPI tasks
    call MPI_Reduce(rhs_intg,rhs_intg_total,5,MPI_REAL,MPI_SUM,master_proc,comm,ierr)

    if (master) then
        print*, '    Integrated RHS rho   = ' , rhs_intg_total(1)
        print*, '    Integrated RHS rho*u = ' , rhs_intg_total(2)
        print*, '    Integrated RHS rho*v = ' , rhs_intg_total(3)
        print*, '    Integrated RHS rho*w = ' , rhs_intg_total(4)
        print*, '    Integrated RHS E     = ' , rhs_intg_total(5)

        if (writetofile) then
            if (step == t1) then
                WRITE(statsUnit,'(6A20)') '#         Time (s)','RHS rho','RHS rho*u','RHS rho*v','RHS rho*w','RHS E'
            end if
            WRITE(statsUnit,'(13ES20.8)') step*dt,rhs_intg_total(1),rhs_intg_total(2),rhs_intg_total(3),rhs_intg_total(4),rhs_intg_total(5)
            CLOSE(statsUnit)
        end if
    end if


    deallocate(     x_flux )
    deallocate(     y_flux )
    deallocate(     z_flux )
    deallocate(        rhs )
    deallocate(       uder )
    deallocate(       vder )
    deallocate(       wder )
    deallocate(       Tder )
    deallocate(         Em )
    deallocate( dilatation )
    deallocate(        tmp )

end subroutine mypostprocess
