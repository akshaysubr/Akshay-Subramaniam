module mpi

    implicit none

    include 'mpif.h'
    
    integer, parameter :: rkind=4                              ! Default kind parameter

    integer :: comm                               ! MPI_COMM_WORLD
    integer, parameter :: stat_size=MPI_STATUS_SIZE
    integer, parameter :: datatype=MPI_REAL
    integer, parameter :: master_proc=0           ! Rank of the master processor

    logical :: master=.FALSE.                     ! Proc is master?
    logical :: x1proc=.FALSE.                     ! Proc is on left x boundary?
    logical :: xnproc=.FALSE.                     ! Proc is on right x boundary?
    logical :: y1proc=.FALSE.                     ! Proc is on left y boundary?
    logical :: ynproc=.FALSE.                     ! Proc is on right y boundary?
    logical :: z1proc=.FALSE.                     ! Proc is on left z boundary?
    logical :: znproc=.FALSE.                     ! Proc is on right z boundary?

    integer :: nprocs                             ! Total number of processors
    integer :: npx                                ! Total no of x processors
    integer :: npy                                ! Total no of y processors
    integer :: npz                                ! Total no of z processors
    integer :: proc                               ! Rank of this processor
    integer :: xproc                              ! x rank of this processor
    integer :: yproc                              ! y rank of this processor
    integer :: zproc                              ! z rank of this processor
    integer :: proc_xr                            ! Proc on right in x-direction
    integer :: proc_xl                            ! Proc on left in x-direction
    integer :: proc_yr                            ! Proc on right in y-direction
    integer :: proc_yl                            ! Proc on left in y-direction
    integer :: proc_zr                            ! Proc on right in z-direction
    integer :: proc_zl                            ! Proc on left in z-direction

    integer :: bcount_x
    integer :: bcount_y
    integer :: bcount_z

    integer :: reqSend_xl,reqSend_xr,reqSend_yl,reqSend_yr,reqSend_zl,reqSend_zr
    integer :: reqRecv_xl,reqRecv_xr,reqRecv_yl,reqRecv_yr,reqRecv_zl,reqRecv_zr

    real(kind=rkind), dimension(:,:,:,:), allocatable :: sendBuf_xl, sendBuf_yl, sendBuf_zl
    real(kind=rkind), dimension(:,:,:,:), allocatable :: sendBuf_xr, sendBuf_yr, sendBuf_zr
    real(kind=rkind), dimension(:,:,:,:), allocatable :: recvBuf_xl, recvBuf_yl, recvBuf_zl
    real(kind=rkind), dimension(:,:,:,:), allocatable :: recvBuf_xr, recvBuf_yr, recvBuf_zr

contains

    function GetProcID(xp,yp,zp) result(p)
        implicit none
        integer, intent(in) :: xp,yp,zp     ! Zero-based indexing as used by MPI
        integer :: p

        p = (zp*npy + yp)*npx + xp

    end function GetProcID

end module mpi
