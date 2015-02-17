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

    use globals, only: u,v,w,rho,p,x_c,y_c,z_c,iodata
    use globals, only: dx,dy,dz,nx,ny,nz,px,py,pz,ax,ay,az
    use globals, only: t1,tf,dt,flen,jobdir
    use interfaces, only: ddx,ddy,ddz,integrate,grad,crossprod
    use functions, only: calc_vorticity
    implicit none
    integer, intent(in) :: step
    integer :: xp,yp,zp
    real(kind=4) :: vortx_intg,vorty_intg,vortz_intg,vortm_intg       ! Integrated rho*omega
    real(kind=4) :: vortx_pres,vorty_pres,vortz_pres,vortm_pres       ! Baroclinic pressure torque
    real(kind=4) :: vortx_visc,vorty_visc,vortz_visc,vortm_visc       ! Viscous torque
    real(kind=4) :: vortx_comp,vorty_comp,vortz_comp,vortm_comp       ! Compressibility effect
    real(kind=4), dimension(:,:,:), allocatable :: dilatation,tmp
    real(kind=4), dimension(:,:,:,:), allocatable :: vort,gradrho,gradp
    character(len=flen) :: statsfile
    integer, parameter :: statsUnit=37
    integer, parameter :: writetofile=0

    if (writetofile == 1) then
        WRITE(statsfile,'(2A)') TRIM(jobdir),'/vortstats.dat'
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
    allocate( dilatation(SIZE(u,1),SIZE(u,2),SIZE(u,3))   )
    allocate(        tmp(SIZE(u,1),SIZE(u,2),SIZE(u,3))   )

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

        ! Calc dilatation
        call ddx(u,dilatation)
        call ddy(v,tmp)
        dilatation = dilatation + tmp
        call ddz(w,tmp)
        dilatation = dilatation + tmp

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
        
        ! ---------------- Vortex gen metrics ----------------
      end do
     end do
    end do

    print*, '    Integrated rho*vorticity = ', vortx_intg, vorty_intg, vortz_intg
    print*, '    Integrated baroclinic pressure term = ', vortx_pres, vorty_pres, vortz_pres
    print*, '    Integrated vortex dilatation term = ', vortx_comp, vorty_comp, vortz_comp

    if (writetofile == 1) then
        if (step == t1) then
            WRITE(statsUnit,'(12A20)') '#         Time (s)','Integrated X vort','Integrated Y vort','Integrated Z vort', &
                                     & 'Max X vort','Max Y vort','Max Z vort','Max vort mag', &
                                     & 'Int X vort gen','Int Y vort gen','Int Z vort gen','Int vort gen mag'
        end if
        ! WRITE(statsUnit,'(12ES20.8)') step*dt,vortx_int,vorty_int,vortz_int, &
        !                          & vortx_max,vorty_max,vortz_max,vortmag_max, &
        !                          & vortx_gen,vorty_gen,vortz_gen,vortmag_gen
        CLOSE(statsUnit)
    end if


    deallocate( vort )
    deallocate( gradrho )
    deallocate( gradp   )
    deallocate( dilatation )
    deallocate( tmp )

end subroutine mypostprocess
