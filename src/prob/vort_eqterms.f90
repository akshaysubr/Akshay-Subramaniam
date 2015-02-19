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

    use globals, only: u,v,w,rho,p,mu,bulk,x_c,y_c,z_c,iodata
    use globals, only: dx,dy,dz,nx,ny,nz,px,py,pz,ax,ay,az
    use globals, only: t1,tf,dt,flen,jobdir
    use interfaces, only: ddx,ddy,ddz,integrate,grad,crossprod
    use functions, only: calc_vorticity
    implicit none
    integer, intent(in) :: step
    integer :: xp,yp,zp,i
    real(kind=4) :: vortx_intg,vorty_intg,vortz_intg,vortm_intg       ! Integrated rho*omega
    real(kind=4) :: vortx_pres,vorty_pres,vortz_pres,vortm_pres       ! Baroclinic pressure torque
    real(kind=4) :: vortx_visc,vorty_visc,vortz_visc,vortm_visc       ! Viscous torque
    real(kind=4) :: vortx_comp,vorty_comp,vortz_comp,vortm_comp       ! Compressibility effect
    real(kind=4), dimension(:,:,:), allocatable :: dilatation,tmp
    real(kind=4), dimension(:,:,:,:), allocatable :: vort,gradrho,gradp,uder,vder,wder
    character(len=flen) :: statsfile
    integer, parameter :: statsUnit=37
    logical, parameter :: writetofile=.TRUE.

    if (writetofile) then
        WRITE(statsfile,'(2A)') TRIM(jobdir),'/vort_eqterms.dat'
        if (step == t1) then
            OPEN(UNIT=statsUnit,FILE=statsfile,FORM='FORMATTED',STATUS='NEW')
        else
            OPEN(UNIT=statsUnit,FILE=statsfile,FORM='FORMATTED',POSITION='APPEND')
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

    ! Proc loop
    do xp=0,px-1
     do yp=0,py-1
      do zp=0,pz-1
        ! Read in grid and data
        call read_thisProc_grid(xp,yp,zp)
        call read_thisProc(step,xp,yp,zp)

        ! ---------------- Vorticity metrics ----------------
        call calc_vorticity(vort)
        vort(:,:,:,1) = rho*vort(:,:,:,1)
        vort(:,:,:,2) = rho*vort(:,:,:,2)
        vort(:,:,:,3) = rho*vort(:,:,:,3)

        ! Integrated rho*omega
        vortx_intg = vortx_intg + integrate(vort(:,:,:,1)) 
        vorty_intg = vorty_intg + integrate(vort(:,:,:,2))
        vortz_intg = vortz_intg + integrate(vort(:,:,:,3))
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

        ! Now get the viscous stress tensor (Put rows in uder, vder, wder)
        do i=1,3
            uder(:,:,:,i) = 2.0D0*mu*uder(:,:,:,i) 
            vder(:,:,:,i) = 2.0D0*mu*vder(:,:,:,i) 
            wder(:,:,:,i) = 2.0D0*mu*wder(:,:,:,i) 
        end do
        uder(:,:,:,1) = uder(:,:,:,1) - (2.0D0*mu/3.0D0 - bulk)*dilatation
        vder(:,:,:,1) = vder(:,:,:,1) - (2.0D0*mu/3.0D0 - bulk)*dilatation
        wder(:,:,:,1) = wder(:,:,:,1) - (2.0D0*mu/3.0D0 - bulk)*dilatation

        ! Compressibility effect
        vortx_comp = vortx_comp - integrate(vort(:,:,:,1)*dilatation)
        vorty_comp = vorty_comp - integrate(vort(:,:,:,2)*dilatation)
        vortz_comp = vortz_comp - integrate(vort(:,:,:,3)*dilatation)

        ! Baroclinic vorticity generation = grad(rho) x grad(p)/(rho)
        call crossprod(vort,gradrho,gradp)

        ! vort now has [ grad(rho) x grad(p) ]
        vortx_pres = vortx_pres + integrate(vort(:,:,:,1)/rho)
        vorty_pres = vorty_pres + integrate(vort(:,:,:,2)/rho)
        vortz_pres = vortz_pres + integrate(vort(:,:,:,3)/rho)

        ! Put div(tau) in gradp
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

        ! Viscous baroclinic torque = grad(rho) x div(tau)/rho
        call crossprod(vort,gradrho,gradp)

        ! vort now has [ grad(rho) x div(tau) ]
        vortx_visc = vortx_visc + integrate(vort(:,:,:,1)/rho)
        vorty_visc = vorty_visc + integrate(vort(:,:,:,2)/rho)
        vortz_visc = vortz_visc + integrate(vort(:,:,:,3)/rho)
        
        ! ---------------- Vortex gen metrics ----------------
      end do
     end do
    end do

    print*, '    Integrated rho*vorticity = ', vortx_intg, vorty_intg, vortz_intg
    print*, '    Integrated baroclinic pressure term = ', vortx_pres, vorty_pres, vortz_pres
    print*, '    Integrated baroclinic viscous term = ', vortx_visc, vorty_visc, vortz_visc
    print*, '    Integrated vortex dilatation term = ', vortx_comp, vorty_comp, vortz_comp

    if (writetofile) then
        if (step == t1) then
            WRITE(statsUnit,'(13A20)') '#         Time (s)','X rho*vort','Y rho*vort','Z rho*vort', &
                                                       & 'X baro pres','Y baro pres','Z baro pres', &
                                                       & 'X baro visc','Y baro visc','Z baro visc', &
                                                       & 'X baro comp','Y baro comp','Z baro comp'
        end if
        WRITE(statsUnit,'(13ES20.8)') step*dt,vortx_intg,vorty_intg,vortz_intg, &
                                            & vortx_pres,vorty_pres,vortz_pres, &
                                            & vortx_visc,vorty_visc,vortz_visc, &
                                            & vortx_comp,vorty_comp,vortz_comp
        CLOSE(statsUnit)
    end if


    deallocate( vort )
    deallocate( gradrho )
    deallocate( gradp )
    deallocate( uder )
    deallocate( vder )
    deallocate( wder )
    deallocate( dilatation )

end subroutine mypostprocess
