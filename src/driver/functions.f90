module functions

contains

    subroutine calc_vorticity(vort)
    
        use globals, only: u,v,w,nx,ny,nz
        use interfaces, only: ddx,ddy,ddz
        implicit none
        real(kind=4), dimension(:,:,:,:) :: vort
        real(kind=4), dimension(:,:,:), allocatable :: tmp,dum

        if(.not. allocated(tmp) ) allocate( tmp(SIZE(vort,1),SIZE(vort,2),SIZE(vort,3)) )
        if(.not. allocated(dum) ) allocate( dum(SIZE(vort,1),SIZE(vort,2),SIZE(vort,3)) )

        call ddy(w,tmp)
        call ddz(v,dum)
        vort(:,:,:,1) = tmp-dum
    
        call ddz(u,tmp)
        call ddx(w,dum)
        vort(:,:,:,2) = tmp-dum
    
        call ddx(v,tmp)
        call ddy(u,dum)
        vort(:,:,:,3) = tmp-dum

        deallocate( tmp )
        deallocate( dum )
    
    end subroutine calc_vorticity

end module
