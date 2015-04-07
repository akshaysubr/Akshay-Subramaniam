module functions

contains

    subroutine calc_vorticity(vort)
    
        use globals, only: rkind
        use globals, only: u,v,w
        use interfaces, only: curl
        implicit none
        real(kind=rkind), dimension(:,:,:,:) :: vort

        call curl(vort(:,:,:,1),vort(:,:,:,2),vort(:,:,:,3),u,v,w)
    
    end subroutine calc_vorticity

end module
