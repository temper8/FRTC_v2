module driver_module
    implicit none

    integer  :: irs
    !common /abcd/ irs
    integer  :: iabsorp
    !common /abcdg/ iabsorp
        
contains

subroutine driver4(ystart,x1,x2,rexi,hmin, derivs)
    use constants
    use runge_kutta_module
    use plasma
    use rt_parameters
    !use trajectory, only: irs, iabsorp
    use dispersion_module
    implicit none
    real(wp), intent(inout)  :: ystart(:)
    real(wp), intent(inout)  :: x1,x2
    real(wp), intent(in)  :: rexi, hmin
    !implicit real*8 (a-h,o-z)
    !procedure (Iderivs_func), pointer, intent(in) :: derivs 
    procedure (Iderivs_func) :: derivs 
    !external derivs
    !common /abcd/ irs
    !common /abcde/ izn!,iw
    !common /abcdg/ iabsorp
    !common /bdeo/ ivar
    !common /bcef/ ynz,ynpopq
    !common /df/ pdec14,pdec24,pdec34,idec
    real(wp) pintld4,pintcl4,pintal4
    common /dg/ pintld4,pintcl4,pintal4
    integer,  parameter :: iturns=1, maxat=3, nvar=4
    real(wp), parameter ::  hbeg=1.d-4 !sav2008
    real(wp)  :: x, xnr, prt, prm, dyd, hnext
    real(wp)  :: yscal(nvar),y(nvar),dydx(nvar),yold(nvar)
    real(wp)  :: eps1, rbord1, hdid, xold, rmm, h
    real(wp)  :: hdrob1, pdec14zz, pdec24zz, pdec34zz
    integer   :: ipr1, iat, i, ii, nstp
    ipr1=0
    iat=0
    x=zero
    eps1=eps
    hdrob1=hdrob
    rbord1=rbord
    hdid=zero
    pintld4=zero
    pintcl4=zero
    pintal4=zero
    pdec14zz=zero
    pdec24zz=zero
    pdec34zz=zero
    xold=x
    do i=1,nvar
        y(i)=ystart(i)
        yold(i)=y(i)
    end do
    rmm=1d+10*irs
    !sav2008
    !old      rexi1=rexi+rrange
    !old      rexi2=rexi-rrange
    !old      if(rexi1.gt.0.95d0) rexi1=1.d10
    !old      if(rexi2.lt.0.05d0) rexi2=-1.d10
    !est      if(rexi1.gt.0.9d0) rexi1=1.1d0
10      continue
    !c--------------------------------------
    !c start integration
    !c--------------------------------------
    do nstp=1,maxstep4
        idec=iturns
        call derivs(x,y,dydx)
        idec=0
        pintld4 = pintld4 + abs((pdec14+pdec14zz)/2d0*hdid)
        pintcl4 = pintcl4 + abs((pdec24+pdec24zz)/2d0*hdid)
        pintal4 = pintal4 + abs((pdec34+pdec34zz)/2d0*hdid)
        pdec14zz=pdec14
        pdec24zz=pdec24
        pdec34zz=pdec34
        if(nstp.eq.1) then
            h=hbeg
            !!var        if(dabs(dydx(3)).ne.zero) h=dabs(hmin1/dydx(3))/hdrob1
            if(dabs(dydx(3)).ne.zero) h=0.5d0*dabs(rrange/dydx(3))/hdrob1
        end if
20          continue
        if(y(3).ge.rbord1.and.dydx(3).gt.zero) then
            !c--------------------------------------
            !c forced reflection from periphery
            !c--------------------------------------
            ivar=3
            izn=-izn
            call disp2(y(3),y(2),y(1),xnr,prt,prm)
            if(ivar.eq.-1) then !out of dispersion curve - restart
                do i=1,nvar
                    y(i)=ystart(i)
                end do
                x=zero
                iat=iat+1
                if(iat.gt.maxat) then
                    if(ipri.gt.1) write (*,*)'turn in driver4 failed'
                    goto 40
                end if
                eps1=eps1/2.d0
                hdrob1=hdrob1*2.d0
                ivar=0
                goto 10
            end if
            irs=-irs
            y(4)=xnr
            call derivs(x,y,dydx)
            if(dydx(3).gt.zero.and.ipri.gt.1) then
                write(*,*)'Unsuccesful turn: r, drds=',y(3),dydx(3)
            end if
            ivar=0
            iat=0
        end if
        !sav2008       if((y(3).gt.rexi1.or.y(3).lt.rexi2)) then  ! exit
        !!    if(dabs(y(3)-rexi).gt.rrange.or.nstp.eq.maxstep4) then  ! exit !sav2008
        if(dabs(y(3)-rexi).gt.rrange) then  ! exit !sav2008
            if(dydx(3).gt.zero) irs=-1
            if(dydx(3).lt.zero) irs=1
            if(dydx(3).eq.zero) then !sav2008
                write(*,*)'exception dr/ds=0 in driver4'
                pause 'zmi na pedal'
                go to 1
            end if
            x2=x
            x1=rmm
            do i=1,nvar
                ystart(i)=y(i)
            end do
            return
        end if
1         continue
        !c---------------------------------------
        !c remember old values
        !c---------------------------------------
        xold=x
        do i=1,nvar
            dyd=dabs(dydx(i))
            yscal(i)=dabs(y(i))+dabs(h*dyd)+1.d-30/(1.d0+dyd)+1.d-30
            yold(i)=y(i)
        end do
        if (y(3)*irs.lt.rmm*irs) rmm=y(3)
30          continue
        
        call runge_kutta_qs(y, dydx, nvar, x, h, eps1, yscal, hdid, hnext, derivs)
        if(y(3).ge.1.d0) then  ! crossed plasma boundary
            do ii=1,nvar
                y(ii)=yold(ii)
            end do
            x=xold
            ipr1=ipr1+1
            if (ipr1.lt.maxat) then
                h=h/3.d0
                goto 30
            end if
            rbord1=y(3)-1.d-4
            goto 20
        end if
        ipr1=0
        if(dabs(hnext).lt.hmin) then
            if(ipri.gt.1) write(*,*)'error in dr4: step is too small'
            goto 40
        end if
        h=hnext
    end do
    if(ipri.gt.1) write(*,*)'error in dr4: too many steps.'
    if(ipri.gt.1) write(*,*)'tet=',y(1),'xm=',y(2),'xend=',y(3)
40    iabsorp=-1
end    
end module driver_module