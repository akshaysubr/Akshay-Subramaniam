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

    use globals, only: rkind
    use globals, only: u,v,w,rho,p,mu,bulk,x_c,y_c,z_c,iodata
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
    real(kind=rkind) :: vortx_intg,vorty_intg,vortz_intg,vortm_intg                               ! Integrated rho*omega
    real(kind=rkind) :: vortx_pres,vorty_pres,vortz_pres,vortm_pres                               ! Baroclinic pressure torque
    real(kind=rkind) :: vortx_visc,vorty_visc,vortz_visc,vortm_visc                               ! Viscous torque
    real(kind=rkind) :: vortx_comp,vorty_comp,vortz_comp,vortm_comp                               ! Compressibility effect
    real(kind=rkind) :: vortx_intg_total,vorty_intg_total,vortz_intg_total,vortm_intg_total       ! Total Integrated rho*omega
    real(kind=rkind) :: vortx_pres_total,vorty_pres_total,vortz_pres_total,vortm_pres_total       ! Total Baroclinic pressure torque
    real(kind=rkind) :: vortx_visc_total,vorty_visc_total,vortz_visc_total,vortm_visc_total       ! Total Viscous torque
    real(kind=rkind) :: vortx_comp_total,vorty_comp_total,vortz_comp_total,vortm_comp_total       ! Total Compressibility effect
    real(kind=rkind), dimension(:,:,:), allocatable :: dilatation,tmp
    real(kind=rkind), dimension(:,:,:,:), allocatable :: vort,gradrho,gradp,uder,vder,wder
    character(len=flen) :: statsfile
    integer, parameter :: statsUnit=37
    logical, parameter :: writetofile=.FALSE.

    if (master) then
        if (writetofile) then
            WRITE(statsfile,'(2A)') TRIM(jobdir),'/vort_eqterms_parallel.dat'
            if (step == t1) then
                OPEN(UNIT=statsUnit,FILE=statsfile,FORM='FORMATTED',STATUS='REPLACE')
            else
                OPEN(UNIT=statsUnit,FILE=statsfile,FORM='FORMATTED',POSITION='APPEND')
            end if
        end if
    end if

    ! Allocate required intermediate arrays
    allocate(       vort(SIZE(u,1),SIZE(u,2),SIZE(u,3),3) )
    allocate(    gradrho(SIZE(u,1),SIZE(u,2),SIZE(u,3),3) )
    allocate(      gradp(SIZE(u,1),SIZE(u,2),SIZE(u,3),3) )
    allocate(       uder(SIZE(u,1),SIZE(u,2),SIZE(u,3),3) )
    allocate(       vder(SIZE(u,1),SIZE(u,2),SIZE(u,3),3) )
    allocate(       wder(SIZE(u,1),SIZE(u,2),SIZE(u,3),3) )
    allocate( dilatation(SIZE(u,1),SIZE(u,2),SIZE(u,3))   )

    vortx_intg = 0.0D0
    vorty_intg = 0.0D0
    vortz_intg = 0.0D0
    vortm_intg = 0.0D0

    vortx_pres = 0.0D0
    vorty_pres = 0.0D0
    vortz_pres = 0.0D0
    vortm_pres = 0.0D0

    vortx_visc = 0.0D0
    vorty_visc = 0.0D0
    vortz_visc = 0.0D0
    vortm_visc = 0.0D0

    vortx_comp = 0.0D0
    vorty_comp = 0.0D0
    vortz_comp = 0.0D0
    vortm_comp = 0.0D0

    ! ---------------- Synthetic fields for debugging --------------
    !VERIFY u = x_c*y_c
    !VERIFY v = -x_c*y_c
    !VERIFY w = 0.0D0
    !VERIFY rho = 1/x_c
    !VERIFY p = 1/y_c
    !VERIFY mu = 1.0D0
    !VERIFY bulk = 1.0D0
    ! ---------------- Synthetic fields for debugging --------------

    !VERIFY if (master) print*, "x, y = ",x_c(ax1+16,ay1+16,az1+16),y_c(ax1+16,ay1+16,az1+16)
    !VERIFY if (master) print*, "u, v = ",u(ax1+16,ay1+16,az1+16), v(ax1+16,ay1+16,az1+16)
    !VERIFY if (master) print*, "u, v = ",x_c(ax1+16,ay1+16,az1+16)*y_c(ax1+16,ay1+16,az1+16),-x_c(ax1+16,ay1+16,az1+16)*y_c(ax1+16,ay1+16,az1+16)
    !VERIFY if (master) print*, "Eu, Ev = ",MAXVAL(ABS(u - x_c*y_c)),MAXVAL(ABS(v - (-x_c*y_c)))

    ! ---------------- Vorticity metrics ----------------
    call calc_vorticity(vort)
    vort(:,:,:,1) = rho*vort(:,:,:,1)
    vort(:,:,:,2) = rho*vort(:,:,:,2)
    vort(:,:,:,3) = rho*vort(:,:,:,3)

    !VERIFY if (master) print*, "x, y = ",x_c(ax1+16,ay1+16,az1+16),y_c(ax1+16,ay1+16,az1+16)
    !VERIFY if (master) print*, "rho*vort = ",vort(ax1+16,ay1+16,az1+16,:)
    !VERIFY if (master) print*, "(rho*vort)_z = ", -( 1.0D0 + y_c(ax1+16,ay1+16,az1+16)/x_c(ax1+16,ay1+16,az1+16) )

    ! Integrated rho*omega
    vortx_intg = integrate(vort(ax1:axn,ay1:ayn,az1:azn,1)) 
    vorty_intg = integrate(vort(ax1:axn,ay1:ayn,az1:azn,2))
    vortz_intg = integrate(vort(ax1:axn,ay1:ayn,az1:azn,3))
    ! ---------------- Vorticity metrics ----------------
    
    ! ---------------- Vortex gen metrics ----------------
    call grad(rho,gradrho(:,:,:,1),gradrho(:,:,:,2),gradrho(:,:,:,3))
    call grad(  p,  gradp(:,:,:,1),  gradp(:,:,:,2),  gradp(:,:,:,3))

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

    ! Now get the Sij tensor in uder, vder and wder
    uder(:,:,:,2) = 0.5D0*(uder(:,:,:,2)+vder(:,:,:,1))
    vder(:,:,:,1) = uder(:,:,:,2)
    uder(:,:,:,3) = 0.5D0*(uder(:,:,:,3)+wder(:,:,:,1))
    wder(:,:,:,1) = uder(:,:,:,3)
    vder(:,:,:,3) = 0.5D0*(vder(:,:,:,3)+wder(:,:,:,2))
    wder(:,:,:,2) = vder(:,:,:,3)

    dilatation = uder(:,:,:,1)+vder(:,:,:,2)+wder(:,:,:,3)
    
    !VERIFY if (master) print*, "x, y = ",x_c(ax1+16,ay1+16,az1+16),y_c(ax1+16,ay1+16,az1+16)
    !VERIFY if (master) print*, "dilatation = ",dilatation(ax1+16,ay1+16,az1+16)
    !VERIFY if (master) print*, "dilatation = ", ( y_c(ax1+16,ay1+16,az1+16) - x_c(ax1+16,ay1+16,az1+16) )

    ! Now get the viscous stress tensor (Put rows in uder, vder, wder)
    do i=1,3
        uder(:,:,:,i) = 2.0D0*mu*uder(:,:,:,i) 
        vder(:,:,:,i) = 2.0D0*mu*vder(:,:,:,i) 
        wder(:,:,:,i) = 2.0D0*mu*wder(:,:,:,i) 
    end do
    uder(:,:,:,1) = uder(:,:,:,1) - (2.0D0*mu/3.0D0 - bulk)*dilatation
    vder(:,:,:,2) = vder(:,:,:,2) - (2.0D0*mu/3.0D0 - bulk)*dilatation
    wder(:,:,:,3) = wder(:,:,:,3) - (2.0D0*mu/3.0D0 - bulk)*dilatation

    !VERIFY if (master) print*, "x, y = ",x_c(ax1+16,ay1+16,az1+16),y_c(ax1+16,ay1+16,az1+16)
    !VERIFY if (master) print*, "vortz_comp = ", -vort(ax1+16,ay1+16,az1+16,3)*dilatation(ax1+16,ay1+16,az1+16)
    !VERIFY if (master) print*, "vortz_comp = ", ( - x_c(ax1+16,ay1+16,az1+16)**2 + y_c(ax1+16,ay1+16,az1+16)**2 )/x_c(ax1+16,ay1+16,az1+16)
    
    ! Compressibility effect
    vortx_comp = -integrate(vort(ax1:axn,ay1:ayn,az1:azn,1)*dilatation(ax1:axn,ay1:ayn,az1:azn))
    vorty_comp = -integrate(vort(ax1:axn,ay1:ayn,az1:azn,2)*dilatation(ax1:axn,ay1:ayn,az1:azn))
    vortz_comp = -integrate(vort(ax1:axn,ay1:ayn,az1:azn,3)*dilatation(ax1:axn,ay1:ayn,az1:azn))

    ! Baroclinic vorticity generation = grad(rho) x grad(p)/(rho)
    call crossprod(vort,gradrho,gradp)

    !VERIFY if (master) print*, "x, y = ",x_c(ax1+16,ay1+16,az1+16),y_c(ax1+16,ay1+16,az1+16)
    !VERIFY if (master) print*, "vortz_pres = ", vort(ax1+16,ay1+16,az1+16,3)/rho(ax1+16,ay1+16,az1+16)
    !VERIFY if (master) print*, "vortz_pres = ", 1.0D0/( x_c(ax1+16,ay1+16,az1+16) * y_c(ax1+16,ay1+16,az1+16)**2 )
    
    ! vort now has [ grad(rho) x grad(p) ]
    vortx_pres = integrate(vort(ax1:axn,ay1:ayn,az1:azn,1)/rho(ax1:axn,ay1:ayn,az1:azn))
    vorty_pres = integrate(vort(ax1:axn,ay1:ayn,az1:azn,2)/rho(ax1:axn,ay1:ayn,az1:azn))
    vortz_pres = integrate(vort(ax1:axn,ay1:ayn,az1:azn,3)/rho(ax1:axn,ay1:ayn,az1:azn))

    ! Put div(tau) in gradp (use dilatation as temporary variable)
    call ddx(uder(:,:,:,1),gradp(:,:,:,1))
    call ddy(uder(:,:,:,2),dilatation)
    gradp(:,:,:,1) = gradp(:,:,:,1) + dilatation
    call ddz(uder(:,:,:,3),dilatation)
    gradp(:,:,:,1) = gradp(:,:,:,1) + dilatation

    call ddx(vder(:,:,:,1),gradp(:,:,:,2))
    call ddy(vder(:,:,:,2),dilatation)
    gradp(:,:,:,2) = gradp(:,:,:,2) + dilatation
    call ddz(vder(:,:,:,3),dilatation)
    gradp(:,:,:,2) = gradp(:,:,:,2) + dilatation

    call ddx(wder(:,:,:,1),gradp(:,:,:,3))
    call ddy(wder(:,:,:,2),dilatation)
    gradp(:,:,:,3) = gradp(:,:,:,3) + dilatation
    call ddz(wder(:,:,:,3),dilatation)
    gradp(:,:,:,3) = gradp(:,:,:,3) + dilatation

    ! Viscous torque = rho * curl( div(tau) / rho )
    ! call crossprod(vort,gradrho,gradp)
    call curl(vort(:,:,:,1),vort(:,:,:,2),vort(:,:,:,3),gradp(:,:,:,1)/rho,gradp(:,:,:,2)/rho,gradp(:,:,:,3)/rho)

    !VERIFY if (master) print*, "x, y = ",x_c(ax1+16,ay1+16,az1+16),y_c(ax1+16,ay1+16,az1+16)
    !VERIFY if (master) print*, "vortz_visc = ", vort(ax1+16,ay1+16,az1+16,3)*rho(ax1+16,ay1+16,az1+16)
    !VERIFY if (master) print*, "vortz_visc = ", 4.0D0/( 3.0D0 * x_c(ax1+16,ay1+16,az1+16) )
    
    ! vort now has [ curl( div(tau) / rho ) ]
    vortx_visc = integrate(vort(ax1:axn,ay1:ayn,az1:azn,1)*rho(ax1:axn,ay1:ayn,az1:azn))
    vorty_visc = integrate(vort(ax1:axn,ay1:ayn,az1:azn,2)*rho(ax1:axn,ay1:ayn,az1:azn))
    vortz_visc = integrate(vort(ax1:axn,ay1:ayn,az1:azn,3)*rho(ax1:axn,ay1:ayn,az1:azn))
    
    !VERIFY stop
    ! ---------------- Vortex gen metrics ----------------

    call MPI_Reduce(vortx_intg,vortx_intg_total,1,MPI_REAL,MPI_SUM,master_proc,comm,ierr)
    call MPI_Reduce(vorty_intg,vorty_intg_total,1,MPI_REAL,MPI_SUM,master_proc,comm,ierr)
    call MPI_Reduce(vortz_intg,vortz_intg_total,1,MPI_REAL,MPI_SUM,master_proc,comm,ierr)

    call MPI_Reduce(vortx_pres,vortx_pres_total,1,MPI_REAL,MPI_SUM,master_proc,comm,ierr)
    call MPI_Reduce(vorty_pres,vorty_pres_total,1,MPI_REAL,MPI_SUM,master_proc,comm,ierr)
    call MPI_Reduce(vortz_pres,vortz_pres_total,1,MPI_REAL,MPI_SUM,master_proc,comm,ierr)

    call MPI_Reduce(vortx_visc,vortx_visc_total,1,MPI_REAL,MPI_SUM,master_proc,comm,ierr)
    call MPI_Reduce(vorty_visc,vorty_visc_total,1,MPI_REAL,MPI_SUM,master_proc,comm,ierr)
    call MPI_Reduce(vortz_visc,vortz_visc_total,1,MPI_REAL,MPI_SUM,master_proc,comm,ierr)

    call MPI_Reduce(vortx_comp,vortx_comp_total,1,MPI_REAL,MPI_SUM,master_proc,comm,ierr)
    call MPI_Reduce(vorty_comp,vorty_comp_total,1,MPI_REAL,MPI_SUM,master_proc,comm,ierr)
    call MPI_Reduce(vortz_comp,vortz_comp_total,1,MPI_REAL,MPI_SUM,master_proc,comm,ierr)

    if (master) then
        print*, '    Integrated rho*vorticity = '           , vortx_intg_total, vorty_intg_total, vortz_intg_total
        print*, '    Integrated baroclinic pressure term = ', vortx_pres_total, vorty_pres_total, vortz_pres_total
        print*, '    Integrated baroclinic viscous term = ' , vortx_visc_total, vorty_visc_total, vortz_visc_total
        print*, '    Integrated vortex dilatation term = '  , vortx_comp_total, vorty_comp_total, vortz_comp_total

        if (writetofile) then
            if (step == t1) then
                WRITE(statsUnit,'(13A20)') '#         Time (s)','X rho*vort','Y rho*vort','Z rho*vort', &
                                                           & 'X baro pres','Y baro pres','Z baro pres', &
                                                           & 'X baro visc','Y baro visc','Z baro visc', &
                                                           & 'X baro comp','Y baro comp','Z baro comp'
            end if
            WRITE(statsUnit,'(13ES20.8)') step*dt,vortx_intg_total,vorty_intg_total,vortz_intg_total, &
                                                & vortx_pres_total,vorty_pres_total,vortz_pres_total, &
                                                & vortx_visc_total,vorty_visc_total,vortz_visc_total, &
                                                & vortx_comp_total,vorty_comp_total,vortz_comp_total
            CLOSE(statsUnit)
        end if
    end if


    deallocate( vort )
    deallocate( gradrho )
    deallocate( gradp )
    deallocate( uder )
    deallocate( vder )
    deallocate( wder )
    deallocate( dilatation )

end subroutine mypostprocess
