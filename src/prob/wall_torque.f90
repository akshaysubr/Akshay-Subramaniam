program postprocess

    use globals, only: lowmem,x_c,y_c,z_c,t1,tf
    implicit none
    integer :: step

    ! Set the low memory flag to true
    lowmem=.TRUE.

    ! Read input file, metadata file and setup globals
    call setup_postprocess

    ! Loop through time steps
    do step=t1,tf
        print*, '> Processing time step ',step
        call mypostprocess(step)
    end do
    
    call cleanup_postprocess

end program postprocess

subroutine mypostprocess(step)

    use globals, only: u,rho,p,x_c,y_c,z_c,iodata
    use globals, only: dx,dy,dz,nx,ny,nz,px,py,pz,ax,ay,az
    use globals, only: t1,tf,dt,flen,jobdir
    use interfaces, only: ddx,ddy,ddz,integrate,grad
    implicit none
    integer, intent(in) :: step
    integer :: xp,yp,zp
    real(kind=4) :: x_cm,y_cm,z_cm,mass,torq,forcey1,forceyf,forcex1,forcexf
    real(kind=4), dimension(:,:,:,:), allocatable :: vort,gradrho,gradp
    character(len=flen) :: statsfile
    integer, parameter :: statsUnit=37

    WRITE(statsfile,'(2A)') TRIM(jobdir),'/wall_torque.dat'
    if (step == t1) then
        OPEN(UNIT=statsUnit,FILE=statsfile,FORM='FORMATTED',STATUS='NEW')
    else
        OPEN(UNIT=statsUnit,FILE=statsfile,FORM='FORMATTED',POSITION='APPEND')
    end if

    x_cm = 0.0D0
    y_cm = 0.0D0
    z_cm = 0.0D0
    mass = 0.0D0
    torq = 0.0D0
    forcey1 = 0.0D0
    forceyf = 0.0D0
    forcex1 = 0.0D0
    forcexf = 0.0D0

    ! Proc loop
    do xp=0,px-1
     do yp=0,py-1
      do zp=0,pz-1
        ! Read in grid and data
        call read_thisProc_grid(xp,yp,zp)
        call read_thisProc(step,xp,yp,zp)

        ! ---------------- Center of mass ----------------
        mass = mass + integrate(rho)
        x_cm = x_cm + integrate(rho*x_c)
        y_cm = y_cm + integrate(rho*y_c)
        z_cm = z_cm + integrate(rho*z_c)
        ! ---------------- Center of mass ----------------
        
        ! ---------------- Boundary pressure Z torque ----------------

        ! Integrate over Y boundaries
        if (yp == 0) then
            torq = torq + SUM( p(:,1,:)*x_c(:,1,:) )*dx*dz
            forcey1 = forcey1 + SUM( p(:,1,:) )*dx*dz
        else if (yp == py-1) then
            torq = torq - SUM(p(:,ay,:)*x_c(:,ay,:))*dx*dz
            forceyf = forceyf + SUM( p(:,ay,:) )*dx*dz
        end if

        ! Integrate over X boundaries
        if (xp == 0) then
            torq = torq + SUM( p(1,:,:)*y_c(1,:,:) )*dy*dz
            forcex1 = forcex1 + SUM( p(1,:,:) )*dy*dz
        else if (xp == px-1) then
            torq = torq - SUM( p(ax,:,:)*y_c(ax,:,:) )*dy*dz
            forcexf = forcexf + SUM( p(ax,:,:) )*dy*dz
        end if
        ! ---------------- Boundary pressure Z torque ----------------
      end do
     end do
    end do

    x_cm = x_cm/mass
    y_cm = y_cm/mass
    z_cm = z_cm/mass

    torq = torq - x_cm*forcey1 + x_cm*forceyf - y_cm*forcex1 + y_cm*forcexf

    print*,'    mass = ',mass
    print*,'    x_cm = ',x_cm
    print*,'    y_cm = ',y_cm
    print*,'    z_cm = ',z_cm
    print*,'    Torque = ',torq

    if (step == t1) then
        WRITE(statsUnit,'(6A20)') '#         Time (s)','mass','X_CM','Y_CM','Z_CM','Z Torque'
    end if
    WRITE(statsUnit,'(6ES20.8)') step*dt,mass,x_cm,y_cm,z_cm,torq
    CLOSE(statsUnit)

end subroutine mypostprocess
