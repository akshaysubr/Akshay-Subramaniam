program postprocess

    use globals, only: x_c,y_c,z_c,t1,tf
    implicit none
    integer :: step

    ! Read input file, metadata file and setup globals
    call setup_postprocess
    call read_allProcs_grid

    ! Loop through time steps
    do step=t1,tf
        print*, '> Processing time step ',step
        call mypostprocess(step)
    end do
    
    call cleanup_postprocess

end program postprocess

subroutine mypostprocess(step)

    use globals, only: u,v,w,rho,p,x_c,y_c,z_c,iodata
    use globals, only: dx,dy,dz,nx,ny,nz,t1,tf,dt,flen,jobdir
    use interfaces, only: ddx,ddy,ddz,integrate,grad,curl
    use functions, only: calc_vorticity
    implicit none
    integer, intent(in) :: step
    real(kind=4) :: vortx_int,vorty_int,vortz_int,vortmag_int
    real(kind=4) :: vortx_max,vorty_max,vortz_max,vortmag_max
    real(kind=4) :: vortx_gen,vorty_gen,vortz_gen,vortmag_gen
    real(kind=4), dimension(:,:,:,:), allocatable :: vort,gradrho,gradp
    real(kind=4), dimension(:,:,:), allocatable :: vx,vy,vz
    character(len=flen) :: statsfile
    integer, parameter :: statsUnit=37
    integer :: ierr

    ! WRITE(statsfile,'(2A)') TRIM(jobdir),'/vortstats.dat'
    ! if (step == t1) then
    !     OPEN(UNIT=statsUnit,FILE=statsfile,FORM='FORMATTED',STATUS='NEW')
    ! else
    !     OPEN(UNIT=statsUnit,FILE=statsfile,FORM='FORMATTED',POSITION='APPEND')
    ! end if

    ! Read data from all procs
    call read_allProcs(step)

    allocate( vort(SIZE(u,1),SIZE(u,2),SIZE(u,3),3), STAT = ierr )
    allocate( vx(SIZE(u,1),SIZE(u,2),SIZE(u,3)) )
    allocate( vy(SIZE(u,1),SIZE(u,2),SIZE(u,3)) )
    allocate( vz(SIZE(u,1),SIZE(u,2),SIZE(u,3)) )

    ! ---------------- Vorticity metrics ----------------
    call calc_vorticity(vort)

    call curl(vx,vy,vz,u,v,w)

    print*, '    Max vx - vort(:,:,:,1) = ', MAXVAL(vx-vort(:,:,:,1))
    print*, '    Max vy - vort(:,:,:,2) = ', MAXVAL(vy-vort(:,:,:,2))
    print*, '    Max vz - vort(:,:,:,3) = ', MAXVAL(vz-vort(:,:,:,3))

    vortx_int   = integrate(vort(:,:,:,1))
    vorty_int   = integrate(vort(:,:,:,2))
    vortz_int   = integrate(vort(:,:,:,3))
    ! vortmag_int = integrate(SQRT( vort(:,:,:,1)**2+vort(:,:,:,2)**2+vort(:,:,:,3)**2 ))
   
    print*, "Calculated integrated vorticity"

    vortx_max   = MAXVAL(ABS(vort(:,:,:,1)))
    vorty_max   = MAXVAL(ABS(vort(:,:,:,2)))
    vortz_max   = MAXVAL(ABS(vort(:,:,:,3)))
    ! vortmag_max = MAXVAL(SQRT( vort(:,:,:,1)**2+vort(:,:,:,2)**2+vort(:,:,:,3)**2 ))
    ! ---------------- Vorticity metrics ----------------
    
    print*,'    Integrated x-vorticity = ',vortx_int
    print*,'    Integrated y-vorticity = ',vorty_int
    print*,'    Integrated z-vorticity = ',vortz_int

    print*,'    Maximum x-vorticity = ',vortx_max
    print*,'    Maximum y-vorticity = ',vorty_max
    print*,'    Maximum z-vorticity = ',vortz_max
    ! print*,'    Maximum vorticity magnitude = ',vortmag_max

    ! allocate( gradrho(SIZE(u,1),SIZE(u,2),SIZE(u,3),3) )
    ! allocate(   gradp(SIZE(u,1),SIZE(u,2),SIZE(u,3),3) )
    ! ! ---------------- Vortex gen metrics ----------------
    ! call grad(rho,gradrho(:,:,:,1),gradrho(:,:,:,2),gradrho(:,:,:,3))
    ! call grad(  p,  gradp(:,:,:,1),  gradp(:,:,:,2),  gradp(:,:,:,3))

    ! ! Baroclinic vorticity generation = grad(rho) x grad(p)/(rho^2)
    ! vort(:,:,:,1) = (gradrho(:,:,:,2)*gradp(:,:,:,3)-gradrho(:,:,:,3)*gradp(:,:,:,2))/(rho**2)
    ! vort(:,:,:,2) = (gradrho(:,:,:,3)*gradp(:,:,:,1)-gradrho(:,:,:,1)*gradp(:,:,:,3))/(rho**2)
    ! vort(:,:,:,3) = (gradrho(:,:,:,1)*gradp(:,:,:,2)-gradrho(:,:,:,2)*gradp(:,:,:,1))/(rho**2)
    ! 
    ! vortx_gen   = integrate(vort(:,:,:,1))
    ! vorty_gen   = integrate(vort(:,:,:,2))
    ! vortz_gen   = integrate(vort(:,:,:,3))
    ! vortmag_gen = integrate(SQRT( vort(:,:,:,1)**2+vort(:,:,:,2)**2+vort(:,:,:,3)**2 ))
    ! 
    ! print*,'    Integrated x-vorticity generation = ',vortx_gen
    ! print*,'    Integrated y-vorticity generation = ',vorty_gen
    ! print*,'    Integrated z-vorticity generation = ',vortz_gen
    ! print*,'    Integrated vorticity generation magnitude = ',vortmag_gen

    ! if (step == 0) then
    !     WRITE(statsUnit,'(12A20)') 'Time (s)','Integrated X vort','Integrated Y vort','Integrated Z vort', &
    !                              & 'Max X vort','Max Y vort','Max Z vort','Max vort mag', &
    !                              & 'Int X vort gen','Int Y vort gen','Int Z vort gen','Int vort gen mag'
    ! end if
    ! WRITE(statsUnit,'(12ES20.8)') step*dt,vortx_int,vorty_int,vortz_int, &
    !                          & vortx_max,vorty_max,vortz_max,vortmag_max, &
    !                          & vortx_gen,vorty_gen,vortz_gen,vortmag_gen
    ! CLOSE(statsUnit)


    ! ---------------- Vortex gen metrics ----------------
    deallocate( vx )
    deallocate( vy )
    deallocate( vz )
    deallocate( vort )
    ! deallocate( gradrho )
    ! deallocate( gradp   )

end subroutine mypostprocess
