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
    integer :: xp,yp,zp,i
    real(kind=4), dimension(:,:,:,:), allocatable :: vort,gradrho,gradp
    real(kind=4), dimension(:), allocatable :: vort_yzint
    character(len=flen) :: statsfile
    integer, parameter :: statsUnit=37

    WRITE(statsfile,'(2A,I0.4)') TRIM(jobdir),'/vortz_yzint.',step
    OPEN(UNIT=statsUnit,FILE=statsfile,FORM='FORMATTED',STATUS='NEW')
    WRITE(statsUnit,'(2A20)') '#         x_c (cm)','YZ Integ. Z vort'

    allocate( vort_yzint(nx) )
    allocate( vort(SIZE(u,1),SIZE(u,2),SIZE(u,3),3) )

    vort_yzint = 0.0D0

    ! Proc loop
    do xp=0,px-1
     do yp=0,py-1
      do zp=0,pz-1
        ! Read in grid and data
        call read_thisProc_grid(xp,yp,zp)
        call read_thisProc(step,xp,yp,zp)

        ! ---------------- Vorticity metrics ----------------
        call calc_vorticity(vort)

        ! Integrate Z vorticity in y and z directions
        vort_yzint(xp*ax+1:(xp+1)*ax) = vort_yzint(xp*ax+1:(xp+1)*ax) + SUM( SUM(vort(:,:,:,3),3),2 )*dy*dz

        ! ---------------- Vorticity metrics ----------------
        
      end do
     end do
     do i=1,ax
         WRITE(statsUnit,'(2ES20.8)') x_c(i,1,1),vort_yzint(xp*ax+i)
     end do
    end do

    CLOSE(statsUnit)


    deallocate( vort_yzint )
    deallocate( vort )

end subroutine mypostprocess
