subroutine ddx_2ndCentral(f,df)

    use globals, only: dx
    implicit none
    real(kind=4), dimension(:,:,:), intent(in) :: f
    real(kind=4), dimension(SIZE(f,1),SIZE(f,2),SIZE(f,3)), intent(out) :: df
    integer :: nx

    nx = SIZE(f,1)

    df(2:nx-1,:,:) = (f(3:nx,:,:) - f(1:nx-2,:,:)) / (2.0*dx)                          ! Interior
    df(     1,:,:) = (-1.5D0*f( 1,:,:) + 2.0D0*f(   2,:,:) - 0.5D0*f(   3,:,:)) / dx   ! Boundary, one sided forward diff
    df(    nx,:,:) = ( 1.5D0*f(nx,:,:) - 2.0D0*f(nx-1,:,:) + 0.5D0*f(nx-2,:,:)) / dx   ! Boundary, one sided backward diff

end subroutine ddx_2ndCentral

subroutine ddy_2ndCentral(f,df)

    use globals, only: dy
    implicit none
    real(kind=4), dimension(:,:,:), intent(in) :: f
    real(kind=4), dimension(SIZE(f,1),SIZE(f,2),SIZE(f,3)), intent(out) :: df
    integer :: ny

    ny = SIZE(f,2)

    df(:,2:ny-1,:) = (f(:,3:ny,:) - f(:,1:ny-2,:)) / (2.0*dy)                          ! Interior
    df(:,     1,:) = (-1.5D0*f(:, 1,:) + 2.0D0*f(:,   2,:) - 0.5D0*f(:,   3,:)) / dy   ! Boundary, one sided forward diff
    df(:,    ny,:) = ( 1.5D0*f(:,ny,:) - 2.0D0*f(:,ny-1,:) + 0.5D0*f(:,ny-2,:)) / dy   ! Boundary, one sided backward diff

end subroutine ddy_2ndCentral

subroutine ddz_2ndCentral(f,df)

    use globals, only: dz
    implicit none
    real(kind=4), dimension(:,:,:), intent(in) :: f
    real(kind=4), dimension(SIZE(f,1),SIZE(f,2),SIZE(f,3)), intent(out) :: df
    integer :: nz

    nz = SIZE(f,3)

    df(:,:,2:nz-1) = (f(:,:,3:nz) - f(:,:,1:nz-2)) / (2.0*dz)                          ! Interior
    df(:,:,     1) = (-1.5D0*f(:,:, 1) + 2.0D0*f(:,:,   2) - 0.5D0*f(:,:,   3)) / dz   ! Boundary, one sided forward diff
    df(:,:,    nz) = ( 1.5D0*f(:,:,nz) - 2.0D0*f(:,:,nz-1) + 0.5D0*f(:,:,nz-2)) / dz   ! Boundary, one sided backward diff

end subroutine ddz_2ndCentral

subroutine d2x_2ndCentral(f,df)

    use globals, only: dx
    implicit none
    real(kind=4), dimension(:,:,:), intent(in) :: f
    real(kind=4), dimension(SIZE(f,1),SIZE(f,2),SIZE(f,3)), intent(out) :: df
    integer :: nx

    nx = SIZE(f,1)

    df(2:nx-1,:,:) = (f(3:nx,:,:) - 2.0D0*f(2:nx-1,:,:) + f(1:nx-2,:,:)) / (dx*dx)              ! Interior
    df( 1,:,:) = ( 2.0D0*f( 1,:,:) - 5.0D0*f(   2,:,:) + 4.0D0*f(   3,:,:) - 1.0D0*f(   4,:,:)) / (dx*dx)   ! Boundary, one sided forward diff
    df(nx,:,:) = (-2.0D0*f(nx,:,:) + 5.0D0*f(nx-1,:,:) - 4.0D0*f(nx-2,:,:) + 1.0D0*f(nx-3,:,:)) / (dx*dx)   ! Boundary, one sided backward diff

end subroutine d2x_2ndCentral

subroutine d2y_2ndCentral(f,df)

    use globals, only: dy
    implicit none
    real(kind=4), dimension(:,:,:), intent(in) :: f
    real(kind=4), dimension(SIZE(f,1),SIZE(f,2),SIZE(f,3)), intent(out) :: df
    integer :: ny

    ny = SIZE(f,2)

    df(:,2:ny-1,:) = (f(:,3:ny,:) - 2.0D0*f(:,2:ny-1,:) + f(:,1:ny-2,:)) / (dy*dy)              ! Interior
    df(:, 1,:) = ( 2.0D0*f(:, 1,:) - 5.0D0*f(:,   2,:) + 4.0D0*f(:,   3,:) - 1.0D0*f(:,   4,:)) / (dy*dy)   ! Boundary, one sided forward diff
    df(:,ny,:) = (-2.0D0*f(:,ny,:) + 5.0D0*f(:,ny-1,:) - 4.0D0*f(:,ny-2,:) + 1.0D0*f(:,ny-3,:)) / (dy*dy)   ! Boundary, one sided backward diff

end subroutine d2y_2ndCentral

subroutine d2z_2ndCentral(f,df)

    use globals, only: dz
    implicit none
    real(kind=4), dimension(:,:,:), intent(in) :: f
    real(kind=4), dimension(SIZE(f,1),SIZE(f,2),SIZE(f,3)), intent(out) :: df
    integer :: nz

    nz = SIZE(f,3)

    df(:,:,2:nz-1) = (f(:,:,3:nz) - 2.0D0*f(:,:,2:nz-1) + f(:,:,1:nz-2)) / (dz*dz)              ! Interior
    df(:,:, 1) = ( 2.0D0*f(:,:, 1) - 5.0D0*f(:,:,   2) + 4.0D0*f(:,:,   3) - 1.0D0*f(:,:,   4)) / (dz*dz)   ! Boundary, one sided forward diff
    df(:,:,nz) = (-2.0D0*f(:,:,nz) + 5.0D0*f(:,:,nz-1) - 4.0D0*f(:,:,nz-2) + 1.0D0*f(:,:,nz-3)) / (dz*dz)   ! Boundary, one sided backward diff

end subroutine d2z_2ndCentral

function integrate_1storder_all(f) result(integral)

    use globals, only: dx,dy,dz
    implicit none
    real(kind=4), dimension(:,:,:), intent(in) :: f
    real(kind=4) :: integral

    integral = SUM(f)*dx*dy*dz

end function integrate_1storder_all

subroutine grad(f,fx,fy,fz)

    use interfaces, only: ddx,ddy,ddz
    implicit none
    real(kind=4), dimension(:,:,:), intent(in) :: f
    real(kind=4), dimension(:,:,:), intent(out) :: fx,fy,fz

    call ddx(f,fx)
    call ddy(f,fy)
    call ddz(f,fz)

end subroutine grad

subroutine curl(fx,fy,fz,u,v,w)

    use interfaces, only: ddx,ddy,ddz
    implicit none
    real(kind=4), dimension(:,:,:), intent(in) :: u,v,w
    real(kind=4), dimension(:,:,:), intent(out) :: fx,fy,fz
    real(kind=4), dimension(:,:,:), allocatable :: tmp,dum

    if(.not. allocated(tmp) ) allocate( tmp(SIZE(u,1),SIZE(u,2),SIZE(u,3)) )
    if(.not. allocated(dum) ) allocate( dum(SIZE(u,1),SIZE(u,2),SIZE(u,3)) )

    print*, '    > In subroutine curl.'

    call ddy(w,tmp)
    call ddz(v,dum)
    fx = tmp-dum

    call ddz(u,tmp)
    call ddx(w,dum)
    fy = tmp-dum

    call ddx(v,tmp)
    call ddy(u,dum)
    fz = tmp-dum

    deallocate( tmp )
    deallocate( dum )

end subroutine curl

subroutine crossprod(rx,ry,rz,ux,uy,uz,vx,vy,vz)

    implicit none
    real(kind=4), dimension(:,:,:), intent(in) :: ux,uy,uz,vx,vy,vz
    real(kind=4), dimension(:,:,:), intent(out) :: rx,ry,rz

    rx = uy*vz - uz*vy
    ry = uz*vx - ux*vz
    rz = ux*vy - uy*vx

end subroutine crossprod
