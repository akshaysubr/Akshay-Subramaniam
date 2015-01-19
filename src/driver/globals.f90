module globals
    implicit none
    integer, parameter :: flen=120
    logical :: verbose=.FALSE.

    integer :: px
    integer :: py
    integer :: pz
    integer :: nprocs

    integer :: nx
    integer :: ny
    integer :: nz
    integer :: ax
    integer :: ay
    integer :: az
    double precision :: dx
    double precision :: dy
    double precision :: dz

    real(kind=4) :: dt
    integer :: t1
    integer :: tf
    integer :: nsteps

    integer, dimension(:,:), allocatable :: procmap
    integer, dimension(:,:,:), allocatable :: invprocmap

    character(len=flen) :: jobdir

    integer :: ndim
    integer :: nvars
    integer :: ns
    integer, parameter :: uInd=1
    integer, parameter :: vInd=2
    integer, parameter :: wInd=3
    integer, parameter :: rhoInd=4
    integer, parameter :: eInd=5
    integer, parameter :: pInd=6
    integer, parameter :: TInd=7
    integer, parameter :: cInd=8
    integer, parameter :: muInd=9
    integer, parameter :: bulkInd=10
    integer, parameter :: ktcInd=11
    integer, parameter :: DiffInd=12
    integer, parameter :: Y1Ind=13

    integer :: ncoords=3
    integer :: xInd
    integer :: yInd
    integer :: zInd

    real(kind=4), dimension(:,:,:,:), allocatable, target :: iodata
    real(kind=4), dimension(:,:,:)  ,  pointer :: u             ! x velocity component
    real(kind=4), dimension(:,:,:)  ,  pointer :: v             ! y velocity component
    real(kind=4), dimension(:,:,:)  ,  pointer :: w             ! z velocity component
    real(kind=4), dimension(:,:,:)  ,  pointer :: rho           ! density
    real(kind=4), dimension(:,:,:)  ,  pointer :: e             ! internal energy
    real(kind=4), dimension(:,:,:)  ,  pointer :: p             ! pressure
    real(kind=4), dimension(:,:,:)  ,  pointer :: T             ! temperature
    real(kind=4), dimension(:,:,:)  ,  pointer :: c             ! speed of sound
    real(kind=4), dimension(:,:,:)  ,  pointer :: mu            ! shear viscosity
    real(kind=4), dimension(:,:,:)  ,  pointer :: bulk          ! bulk viscosity
    real(kind=4), dimension(:,:,:)  ,  pointer :: ktc           ! thermal conductivity
    real(kind=4), dimension(:,:,:)  ,  pointer :: Diff          ! diffusivity
    real(kind=4), dimension(:,:,:,:),  pointer :: Y             ! mass fractions
    real(kind=4), dimension(:,:,:)  ,  pointer :: x_c           ! x coordinate
    real(kind=4), dimension(:,:,:)  ,  pointer :: y_c           ! y coordinate
    real(kind=4), dimension(:,:,:)  ,  pointer :: z_c           ! z coordinate

    double precision, parameter :: pi=3.14159265358979323846d0,zero=0.0d0,one=1.0d0,two=2.0d0

end module globals

subroutine setup_pointers

    use globals, only: u,v,w,rho,e,p,T,c,mu,bulk,ktc,Diff,Y,x_c,y_c,z_c,iodata,ns
    use globals, only: uInd,vInd,wInd,rhoInd,eInd,pInd,TInd,cInd,muInd,bulkInd,ktcInd,DiffInd,Y1Ind,xInd,yInd,zInd
    implicit none

    u    => iodata(:,:,:,    uInd)
    v    => iodata(:,:,:,    vInd)
    w    => iodata(:,:,:,    wInd)
    rho  => iodata(:,:,:,  rhoInd)
    e    => iodata(:,:,:,    eInd)
    p    => iodata(:,:,:,    pInd)
    T    => iodata(:,:,:,    TInd)
    c    => iodata(:,:,:,    cInd)
    mu   => iodata(:,:,:,   muInd)
    bulk => iodata(:,:,:, bulkInd)
    ktc  => iodata(:,:,:,  ktcInd)
    Diff => iodata(:,:,:, DiffInd)
    Y    => iodata(:,:,:, Y1Ind:Y1Ind+ns)
    x_c  => iodata(:,:,:,    xInd)
    y_c  => iodata(:,:,:,    yInd)
    z_c  => iodata(:,:,:,    zInd)

end subroutine setup_pointers
