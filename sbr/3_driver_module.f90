module driver_module
    use kind_module
    implicit none
    integer, parameter :: length = 5000000

    real(wp) dland(length),dcoll(length),perpn(length),dalf(length)
    real(wp) vel(length),tetai(length)
    real(wp) xnpar(length)
    integer izz(length),iww(length),jrad(length)
    !!common/agh/xnpar,vel,dland,dcoll,dalf,perpn,tetai,jrad,iww,izz
        
    integer  :: irs
    !common /abcd/ irs
    integer  :: iabsorp
    !common /abcdg/ iabsorp

    real(wp) :: rzz,tetzz,xmzz
    !!common /abc/ rzz,tetzz,xmzz    
    integer  :: iznzz,iwzz,irszz
    !!common /abc/ iznzz,iwzz,irszz

    real(wp) :: hrad
    !common /bcg/ hrad    


    integer  :: im4
    !common /bg/ im4    
contains

    subroutine driver2(ystart,x1,x2,xsav,hmin,h1, pabs) !sav2008
        use constants
        use plasma
        use rt_parameters
        !use manager_mod
        !use trajectory
        use dispersion_module
        !use driver_module
        implicit real*8 (a-h,o-z)
        !external extd2
        real*8 pabs
        !common /abc/ rzz,tetzz,xmzz,iznzz,iwzz,irszz
        !common /abcd/ irs
        !common /abcde/ izn!,iw
        !common /abcdg/ iabsorp
        !common /bcef/ ynz,ynpopq
        !common /bcg/ hrad
        !common /cefn/ iconv,irefl
        !common /ceg/ ipow,jfoundr
        integer :: ind
        common /cmn/ ind
        dimension ystart(2)
        integer, parameter :: nvar=2
        dimension yscal(nvar),y(nvar),dydx(nvar),yold(nvar),dyold(nvar)
        integer :: i, ii, irep, nstp
        x=x1
        h=dsign(h1,x2-x1)
        ind=0
        ipow=-1
        xold=x
        hsav=hrad*irs
        hdid=zero
        do i=1,nvar
            y(i)=ystart(i)
            yold(i)=y(i)
        end do
        !c-----------------------------
        !c            start moving
        !c-----------------------------
        do nstp=1,maxstep2
            !c---------------------------------------
            !c netpoint control
            !c---------------------------------------
            dstsav=dabs(x-xsav)
            if(dstsav.lt.tin) then
                ipow=ipow+2
                jfoundr=idnint(x/hrad)
                if(jfoundr.le.0) jfoundr=1
                if(jfoundr.gt.nr) jfoundr=nr
            end if
            call extd2(x,y,dydx)
            irep=0
            if(iconv+irefl.ne.0) then
                !!       if(iconv+irefl.ne.0.or.ynz.lt.0.d0) then
                !c---------------------------------------------
                !c made step to nontransparent zone-return back
                !c----------------------------------------------
                x=xold
                do ii=1,nvar
                    y(ii)=yold(ii)
                    dydx(ii)=dyold(ii)
                end do
                irep=1
                h=hdid/2
                hdid=h
                ipow=0
                if(dabs(h).lt.hmin1) then
                    ind=3
                    go to 20
                end if
                ynz=ynz0
                go to 10
            end if            
            !c--------------------------------------
            !c memorize step data
            !c--------------------------------------
            xold=x
            do i=1,nvar
                dyd=dabs(dydx(i))
                yscal(i)=dabs(y(i))+dabs(h*dyd)+1.d-30/(1.d0+dyd)
                yold(i)=y(i)
                dyold(i)=dydx(i)
            end do
            ynz0=ynz
            if(ipow.gt.0) then !integrate power equation
                call dql1(pabs)
                if(iabsorp.eq.1) then !absorption
                    rzz=x
                    tetzz=y(1)
                    xmzz=y(2)
                    iznzz=izn
                    iwzz=iw
                    irszz=irs
                    return
                end if
                if(iabsorp.eq.-1) return !problem
                ipow=0
                xsav=xsav-hsav
            end if
            !c--------------------------------------
            !c choose step size
            !c--------------------------------------
            dst3=(x-xsav)*(x+h-xsav)
            if(dst3.lt.zero.and.irep.eq.0) h=xsav-x
            if(x.gt.rbord.and.h.gt.zero) then
                ind=2
                go to 20
            end if
10              dst1=(x-rbord)*(x+h-rbord)
            dst2=x*(x+h)
            if((dst1.lt.zero.and.irs.eq.-1).or.dst2.lt.zero) then
                h=h/2.d0
                if(dabs(h).lt.hmin1) then
                    ind=4
                    go to 20
                end if
                go to 10
            end if
            !c--------------------------------------
            !c find solution at x=x+hdid
            !c---------------------------------------
            ynz0=ynz
            call difeq(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,extd2)
20              continue
            if(ind.ne.0) then !exit
                xsav=xsav+hsav
                x2=x
                do i=1,nvar
                    ystart(i)=y(i)
                end do
                ynz=ynz0
                return
            end if
            !c---------------------------------------
            if(dabs(hnext).lt.hmin) then
                if(ipri.gt.1) write(*,*) 'exit driver2: step is too small'
                go to 40
            end if
            h=hnext
        end do
            !c---------------------------------------
        if (ipri.gt.1) write (*,*) 'error in driver2: too many steps'
40      iabsorp=-1
        return
1001    format (10(e14.7,1x))

    end

    subroutine extd2(x,y,dydx)
        use dispersion_module, only: disp2
        implicit real*8 (a-h,o-z)
        dimension y(*),dydx(*)
        tt=y(1)
        xm=y(2)
        call disp2(x,xm,tt,xnr,prt,prm)
        dydx(1)=-prm
        dydx(2)=prt
    end

    subroutine dql1(pabs) !sav2008
        use constants, only: clt, zero
        use rt_parameters
        use plasma, only : fvt, vperp
        use current
        !use spectrum1D, only: pabs
        !use trajectory
        use dispersion_module
        use manager_mod, only: pow, inak, lenstor, lfree
        !use driver_module !, only:  iabsorp
        implicit real*8 (a-h,o-z)
        real*8 pabs
        real*8 radth
        dimension an1(length),an2(length)
        common /xn1xn2/ an1,an2
        common /vth/ vthc(length),poloidn(length)
        common /a0ghp/ vlf,vrt,dflf,dfrt
        !common /abcdg/ iabsorp
        !common /acg/ pow
        !common /ag/ inak,lenstor,lfree
        !common /bcg/ hrad
        !common /bg/ im4
        !common /ceg/ ipow,jfoundr
        !common /eg1/ vfound,ifound
        !common /eg2/ pdec1,pdec2,pdec3,pdecv,pdecal,dfdv,icf1,icf2
        !common /eg3/ cf1,cf2,cf3,cf4,cf5,cf6
        common /dg/ pintld4,pintcl4,pintal4
        integer :: i, j, npoloid, ifast, idir
        powpr=pow
        iabsorp=0
        hdis=hrad
        vz=vfound
        i=ifound
        if(i.eq.0) i=1
        j=jfoundr
        refr=cf1
        tet_i=cf2
        npoloid=cf6
        xparn=cf3
        xan1=cf4
        xan2=cf5
        ifast=icf1
        idir=icf2
        dek3=zero
        dfsr=(vlf*dflf+vrt*dfrt)/2d0*(vrt-vlf)
        vsr=(vrt+vlf)*(vrt-vlf)/2d0
        !c--------------------------------------
        !c   find power
        !c--------------------------------------
        if(im4.eq.1) then
            !!       pintld=-pintld4*dfdv
            pintld=dabs(pintld4*dfdv)
            pintcl=dabs(pintcl4)
            if(itend0.gt.0) then
                argum=clt/(refr*valfa)
                dek3=zatukh(argum,j,vperp,kv)
            end if
            pintal=dabs(pintal4*dek3)
            dcv=pintld4/vsr
        else
            pintld=dabs(pdec1*hdis)
            pintcl=dabs(pdec2*hdis)
            pintal=dabs(pdec3*hdis)
            dcv=pdecv*hdis/vsr
        end if
        if(pabs.ne.zero) then
            powd=pow*dexp(-2d0*pintld)
            powccc=dexp(-2d0*pintcl)
            powcol=powd*powccc
            powal=powcol*dexp(-2d0*pintal)
            pow=powal
        end if
        if(pow.le.pabs) iabsorp=1
        pil=pintld
        pic=pintcl
        pia=pintal
        call dfind(j,i,vz,powpr,pil,pic,pia,dfsr,dcv &
                                ,refr,vlf,vrt,ifast)
        !c-----------------------------------
        !c      memorize trajectory
        !c----------------------------------
        inak=inak+1
        if(inak.eq.lenstor) then
            write(*,*)'storage capacity exceeded !'
            iabsorp=-1
            inak=lenstor-1
            return
        end if
        vel(inak)=vz
        perpn(inak)=refr
        poloidn(inak)=npoloid
        tetai(inak)=tet_i
        radth=dble(j)/dble(31)
        vthc(inak)=3.d10/fvt(radth)
        iww(inak)=ifast
        izz(inak)=idir
        xnpar(inak)=xparn
        an1(inak)=xan1
        an2(inak)=xan2
        if(im4.eq.1) then
            jrad(inak)=-j
            dland(inak)=pintld4
            dcoll(inak)=pintcl4
            dalf(inak)=pintal4
            im4=0
            return
        end if
        jrad(inak)=j
        dland(inak)=pdecv
        dalf(inak)=pdecal
        if(ipow.ne.1) dcoll(inak)=powccc
        if(ipow.eq.1) dcoll(inak)=1d0
    end

!------------------------------------------------------------------------------------------------
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