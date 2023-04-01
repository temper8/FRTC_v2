module iterator_mod
    use kind_module   
    implicit none
    real(wp) :: vmid(100),vz1(100),vz2(100)
    integer  :: ibeg(100),iend(100)

    real(wp) :: vrj(101),dj(101),djnew(1001)
    real(wp) :: dj2(101),d2j(101)

    real(wp), dimension(:), allocatable:: vvj, vdfj

    real(wp) :: vgrid(101,100), dfundv(101,100)
    !!common/gridv/vgrid(101,100),dfundv(101,100)
    integer  :: nvpt
    !!common/gridv/nvpt
    integer :: ipt1, ipt2, ipt
    integer, parameter :: kpt1=20, kpt3=20

    integer :: iterat
    real(wp) :: psum4
    !!common /vvv2/ psum4
    real(wp) plost,pnab
    !!common /a0a4/ plost,pnab
contains
    subroutine init_iteration
        use constants, only : zero
        use rt_parameters, only : itend0
        use current
        use plasma, only: cltn
        implicit none
        ppv1=zero
        ppv2=zero
        psum4=zero
        pnab=zero
        plost=zero
        dql=zero
        dq1=zero
        dq2=zero
        dncount=zero
        vzmin=cltn
        vzmax=-cltn
        pdl=zero
        pdc=zero
        pda=zero
        pdfast=zero
        if(itend0.gt.0) then
              dqi0=zero
        end if
    end 

end module iterator_mod