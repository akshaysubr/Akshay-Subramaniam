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
