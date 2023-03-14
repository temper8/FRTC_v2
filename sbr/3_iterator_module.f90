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
    subroutine recalculate_f_for_a_new_mesh(ispectr)
        !!   recalculate f' for a new mesh
        use constants, only : zero
        use rt_parameters, only : nr, ni1,ni2
        use plasma, only: vt0, fvt, cltn
        use current, only: vzmin, vzmax
        use maxwell, only: i0, vij, dfij
        implicit none
        integer, intent(in) :: ispectr
        
        integer i, j, k
        real(wp) :: cdel, dfout

        integer :: klo,khi,ierr
        real(wp) :: r, hr, vt, vto, vmax
        real(wp) :: v1, v2, vp1, vp2
        hr = 1.d0/dble(nr+1)
        k=(3-ispectr)/2
        do j=1,nr
            r=hr*dble(j)
            vt=fvt(r)
            vto=vt/vt0
            if(iterat.gt.0) then
                v1=dmin1(vzmin(j),vz1(j))
                v2=dmax1(vzmax(j),vz2(j))
            else
                v1=vzmin(j)
                v2=vzmax(j)
            end if
            vmax=cltn/vto
            vp1=v1/vto
            vp2=v2/vto
            call gridvel(vp1,vp2,vmax,cdel,ni1,ni2,ipt1,kpt3,vrj)
            do i=1,i0
                vvj(i)=vij(i,j)
                vdfj(i)=dfij(i,j,k) !=dfundv(i,j)*vto**2
            end do
            do i=1,ipt
                call lock(vvj,i0,vrj(i),klo,khi,ierr)
                if(ierr.eq.1) then
            !!!         if(vrj(i).gt.vvj(i0)) exit
                    write(*,*)'lock error in new v-mesh'
                    write(*,*)'j=',j,' i0=',i0
                    write(*,*)'vvj(1)=',vvj(1),' vvj(i0)=',vvj(i0)
                    write(*,*)'i=',i,' vrj(i)=',vrj(i)
                    write(*,*)
                    pause'next key = stop'
                    stop
                end if
                call linf(vvj,vdfj,vrj(i),dfout,klo,khi)
                vgrid(i,j)=vrj(i)*vto
                dfundv(i,j)=dfout/vto**2
                if(dfundv(i,j).gt.zero) dfundv(i,j)=zero
            end do
            vz1(j)=v1
            vz2(j)=v2
        end do
    end subroutine
end module iterator_mod