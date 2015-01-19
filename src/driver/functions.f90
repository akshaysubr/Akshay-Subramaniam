subroutine calc_vorticity(vort)

    use globals, only: u,v,w
    use interfaces, only: ddx,ddy,ddz
    implicit none
    real(kind=4), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3),3), intent(out) :: vort
    real(kind=4), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: tmp,dum

    call ddy(w,tmp)
    call ddz(v,dum)
    vort(:,:,:,1) = tmp-dum

    call ddz(u,tmp)
    call ddx(w,dum)
    vort(:,:,:,2) = tmp-dum

    call ddx(v,tmp)
    call ddy(u,dum)
    vort(:,:,:,3) = tmp-dum

end subroutine calc_vorticity

