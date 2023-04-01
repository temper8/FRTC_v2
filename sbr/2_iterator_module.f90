module iterator_mod
    use kind_module   
    implicit none
    real(wp) :: vmid(100),vz1(100),vz2(100)
    integer  :: ibeg(100),iend(100)

    real(wp) :: vrj(101),dj(101),djnew(1001)
    real(wp) :: dj2(101),d2j(101)

    real(wp) :: vgrid(101,100), dfundv(101,100)
    !!common/gridv/vgrid(101,100),dfundv(101,100)
    integer  :: nvpt
    !!common/gridv/nvpt
    integer :: ipt1, ipt2, ipt


    integer  :: iterat
    real(wp) :: psum4
    !!common /vvv2/ psum4
    real(wp) ::plost,pnab
    !!common /a0a4/ plost,pnab
contains

end module iterator_mod