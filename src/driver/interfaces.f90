module interfaces

    interface ddx
        subroutine ddx_2ndCentral(f,df)
            use globals, only: dx,rkind
            real(kind=rkind), dimension(:,:,:), intent(in) :: f
            real(kind=rkind), dimension(SIZE(f,1),SIZE(f,2),SIZE(f,3)), intent(out) :: df
        end subroutine ddx_2ndCentral
    end interface

    interface ddy
        subroutine ddy_2ndCentral(f,df)
            use globals, only: dy,rkind
            real(kind=rkind), dimension(:,:,:), intent(in) :: f
            real(kind=rkind), dimension(SIZE(f,1),SIZE(f,2),SIZE(f,3)), intent(out) :: df
        end subroutine ddy_2ndCentral
    end interface

    interface ddz
        subroutine ddz_2ndCentral(f,df)
            use globals, only: dz,rkind
            real(kind=rkind), dimension(:,:,:), intent(in) :: f
            real(kind=rkind), dimension(SIZE(f,1),SIZE(f,2),SIZE(f,3)), intent(out) :: df
        end subroutine ddz_2ndCentral
    end interface

    interface d2x
        subroutine d2x_2ndCentral(f,df)
            use globals, only: dx,rkind
            implicit none
            real(kind=rkind), dimension(:,:,:), intent(in) :: f
            real(kind=rkind), dimension(SIZE(f,1),SIZE(f,2),SIZE(f,3)), intent(out) :: df
        end subroutine d2x_2ndCentral
    end interface

    interface d2y
        subroutine d2y_2ndCentral(f,df)
            use globals, only: dy,rkind
            implicit none
            real(kind=rkind), dimension(:,:,:), intent(in) :: f
            real(kind=rkind), dimension(SIZE(f,1),SIZE(f,2),SIZE(f,3)), intent(out) :: df
        end subroutine d2y_2ndCentral
    end interface

    interface d2z
        subroutine d2z_2ndCentral(f,df)
            use globals, only: dz,rkind
            implicit none
            real(kind=rkind), dimension(:,:,:), intent(in) :: f
            real(kind=rkind), dimension(SIZE(f,1),SIZE(f,2),SIZE(f,3)), intent(out) :: df
        end subroutine d2z_2ndCentral
    end interface

    interface integrate
        function integrate_1storder_all(f) result(integral)
            use globals, only: dx,dy,dz,rkind
            implicit none
            real(kind=rkind), dimension(:,:,:), intent(in) :: f
            real(kind=rkind) :: integral
        end function integrate_1storder_all
    end interface

    interface grad
        subroutine grad(f,fx,fy,fz)
            use globals, only: rkind
            implicit none
            real(kind=rkind), dimension(:,:,:), intent(in) :: f
            real(kind=rkind), dimension(:,:,:), intent(out) :: fx,fy,fz
        end subroutine grad
    end interface

    interface curl
        subroutine curl(fx,fy,fz,u,v,w)
            use globals, only: rkind
            implicit none
            real(kind=rkind), dimension(:,:,:), intent(in) :: u,v,w
            real(kind=rkind), dimension(:,:,:), intent(out) :: fx,fy,fz
        end subroutine curl
    end interface


    interface crossprod
        subroutine crossprod_components(rx,ry,rz,ux,uy,uz,vx,vy,vz)
            use globals, only: rkind
            implicit none
            real(kind=rkind), dimension(:,:,:), intent(in) :: ux,uy,uz,vx,vy,vz
            real(kind=rkind), dimension(:,:,:), intent(out) :: rx,ry,rz
        end subroutine crossprod_components
        subroutine crossprod_arrays(r,u,v)
            use globals, only: rkind
            implicit none
            real(kind=rkind), dimension(:,:,:,:), intent(in) :: u,v
            real(kind=rkind), dimension(:,:,:,:), intent(out) :: r
        end subroutine crossprod_arrays
    end interface

end module interfaces
