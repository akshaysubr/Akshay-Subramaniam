module interfaces

    interface ddx
        subroutine ddx_2ndCentral(f,df)
            use globals, only: dx
            real(kind=4), dimension(:,:,:), intent(in) :: f
            real(kind=4), dimension(SIZE(f,1),SIZE(f,2),SIZE(f,3)), intent(out) :: df
        end subroutine ddx_2ndCentral
    end interface

    interface ddy
        subroutine ddy_2ndCentral(f,df)
            use globals, only: dy
            real(kind=4), dimension(:,:,:), intent(in) :: f
            real(kind=4), dimension(SIZE(f,1),SIZE(f,2),SIZE(f,3)), intent(out) :: df
        end subroutine ddy_2ndCentral
    end interface

    interface ddz
        subroutine ddz_2ndCentral(f,df)
            use globals, only: dz
            real(kind=4), dimension(:,:,:), intent(in) :: f
            real(kind=4), dimension(SIZE(f,1),SIZE(f,2),SIZE(f,3)), intent(out) :: df
        end subroutine ddz_2ndCentral
    end interface

    interface integrate
        function integrate_1storder_all(f) result(integral)
            use globals, only: dx,dy,dz
            implicit none
            real(kind=4), dimension(:,:,:), intent(in) :: f
            real(kind=4) :: integral
        end function integrate_1storder_all
    end interface

    interface grad
        subroutine grad(f,fx,fy,fz)
            implicit none
            real(kind=4), dimension(:,:,:), intent(in) :: f
            real(kind=4), dimension(:,:,:), intent(out) :: fx,fy,fz
        end subroutine grad
    end interface
end module interfaces
