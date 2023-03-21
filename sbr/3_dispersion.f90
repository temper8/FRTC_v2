module dispersion_module
    use kind_module    
    implicit none

    real(wp) :: yn3
    !! common /abefo/ yn3

    integer :: ivar
    !!common /bdeo/ ivar   

    integer :: icall1, icall2
    !common /aef2/ icall1,icall2        

    integer  :: iroot
    !!common /beo/ iroot
    integer  :: izn
    !common /abcde/ izn
    integer  :: ider
    !common /be2/ ider

    real(wp) :: xnr1,xnr2,xnr3,xnr4
    !common /be1/ xnr1,xnr2,xnr3,xnr4

    real(wp) ::ynz, ynpopq
    !common /bcef/ ynz,ynpopq

contains
    subroutine disp2(pa,yn2,ptet,xnro,prt,prm)
        use constants
        use approximation
        use plasma
        use rt_parameters
        !use manager_mod, only: ivar, yn3, icall1, icall2
        !use trajectory !, only: iroot, izn, ynz,ynpopq
        implicit real*8 (a-h,o-z)
        !common /abcde/ izn!,iw
        !common /bcef/ ynz,ynpopq
        !common /aef2/ icall1,icall2
        !common /be1/ xnr1,xnr2,xnr3,xnr4
        !common /be2/ ider
        integer iconv, irefl
        common /cefn/ iconv,irefl

        integer ipow,jfoundr
        common /ceg/ ipow,jfoundr
        
        integer ifound
        common /eg1/ vfound,ifound
        common /eg2/ pdec1,pdec2,pdec3,pdecv,pdecal,dfdv,icf1,icf2
        common /eg3/ cf1,cf2,cf3,cf4,cf5,cf6
        common/fj/dhdm,dhdnr,dhdtet,dhdr,ddn,dhdn3,dhdv2v,dhdu2u
        common/direct/znakstart
        common/metrika/g11,g12,g22,g33,gg,g,si,co
        common/fjham/ham
        integer  jr, icf1, icf2

        iconv=0
        irefl=0
        if(pa.ge.one.or.pa.le.zero) goto 70
        icall1=icall1+1
        xdl=fdf(pa,cdl,ncoef,xdlp)
        xly=fdf(pa,cly,ncoef,xlyp)
        xgm=fdf(pa,cgm,ncoef,xgmp)
        xmy=fdf(pa,cmy,ncoef,xmyp)
        xlyv=xlyp*pa+xly
        cotet=dcos(ptet)
        sitet=dsin(ptet)
        dxdr=-xdlp+cotet-xgmp*sitet**2
        dxdt=-(pa+two*xgm*cotet)*sitet
        dzdr=xlyv*sitet
        dzdt=xly*pa*cotet
        x0=r0/rm-xdl+pa*cotet-xgm*sitet**2
        dxdrdt=-sitet-two*xgmp*sitet*cotet
        dzdrdt=xlyv*cotet
        dxdtdt=-pa*cotet-two*xgm*(cotet**2-sitet**2)
        dzdtdt=-xly*pa*sitet
        x0t=dxdt
        !--------------------------------------
        ! components of metric tensor
        !--------------------------------------
        g11=dxdr**2+dzdr**2
        g22=dxdt**2+dzdt**2
        g12=dxdr*dxdt+dzdr*dzdt
        g33=x0**2
        xj=(dzdr*dxdt-dxdr*dzdt)**2  !gg=g11*g22-g12*g12
        gg=xj
        g=xj*g33
        g2v1=one/dsqrt(g22)
        g2jq=dsqrt(g22/xj)
        g3v=one/dsqrt(g33)
        !--------------------------------------
        !  magnetic field
        !--------------------------------------
        bt=b_tor*(r0/rm)/x0
        bp=g2jq*g3v*xmy
        b=dsqrt(bp*bp+bt*bt)
        si=bp/b
        co=bt/b
        if(ivar.eq.1) return
        !---------------------------------------
        ! components of dielectric tensor
        !---------------------------------------
        !sav2008      pn=fn(pa)
        !var      pn=fn1(pa,fnr)
        !!      pn=fn1(pa,fnr)
        !!      pn=fn2(pa,fnr,fnrr) !sav2008
        if(inew.eq.0) then !vardens
            pn=fn1(pa,fnr)
        else
            pn=fn2(pa,fnr,fnrr)
        end if
        wpq=c0**2*pn
        whe=b*c1
        v=wpq/ww**2
        u1=whe/ww
        u=u1**2
        e1=one-v*(one/xmi-one/u)
        e2=v/u1
        e3=one-v
        !-------------------------------------
        ! dispersion equation
        !--------------------------------------
        !sav2008      if(ivar.eq.2) yn2=(ynz-yn3*co*g3v)/(si*g2v1)
        ynz=yn2*si*g2v1+yn3*co*g3v
        ynzq=ynz**2
        as=e1
        bs=-(e1**2-e2**2+e1*e3-(e1+e3)*ynzq)
        cs=e3*(e1**2-e2**2-two*e1*ynzq+ynzq**2)
        !-----------------------------
        !est !sav2009
        pnew=zero
        yny= - (yn2*g2v1*co-yn3*g3v*si)
        if(inew.gt.0) then
            if(inew.eq.1) then
                    yny= - (yn2*g2v1*co-yn3*g3v*si)
                else if(inew.eq.2) then
                    yny= - g2jq*(yn2*g2v1*co-yn3*g3v*si)
            end if
            gpr=c0**2/ww**2/u1*fnr*xsz
            pnew=yny*gpr
            bs=bs+pnew
            cs=cs+pnew*(ynzq-e3)
        end if
        !------------------------------------
        dls=bs*bs-4d0*as*cs

        !c      write(*,*)'rho=',pa,' teta=',ptet
        !c      write(*,*)'N2=',yn2,' N3=',yn3
        !c      write(*,*)'Npar=',ynz,' e1=',e1
        !c      write(*,*)'v=',v,' u=',u
        !c      write(*,*)'whe=',whe,' ww=',ww
        !c      write(*,*)'e2=',e2,' e3=',e3
        !c      write(*,*)'bs=',bs,' as=',as
        !c      write(*,*)'cs=',cs,' dls=',dls
        !c      pause
        
        if(dls.lt.zero) then
            goto (60,20,10) iroot
10          xnr1=1d+10
            xnr2=1d+10
            xnr3=1d+10
            xnr4=1d+10
            return
20          prt=dls
            prm=666d0
            return
        end if
30      continue
        dl1=dfloat(iw)*dsqrt(dls)/two/as
        if(iw.eq.-1) ynpopq=-bs/(two*as)+dl1
        if(iw.eq.1)  ynpopq=two*cs/(-bs-two*as*dl1)
        if(iroot.eq.3) ynpopq1=-bs/(two*as)-dl1

        !cc      write(*,*)'iw=',iw,' izn=',izn,' Nperp=',dsqrt(ynpopq)
        !cc      write(*,*)'Nperp2=',ynpopq,' ynpopq1=',-bs/(two*as)-dl1
        !cc      pause

        if(ynpopq.lt.zero.and.iroot.eq.1) goto 70
        al=g22/xj
        bl=-yn2*g12/xj
        cl=g11*yn2**2/xj+yn3**2/g33-ynzq-ynpopq
        if(iroot.eq.3) cl1=g11*yn2**2/xj+yn3**2/g33-ynzq-ynpopq1
        dll=bl*bl-al*cl
        if(iroot.eq.2) then
            prt=dls
            prm=dll
            if(dll.ge.zero) then !sav2008
                !!old variant:
                !cc        dl2=-dfloat(izn)*dsqrt(dll)/al
                !cc        if(izn.eq.1) xnr=-bl/al+dl2
                !cc        if(izn.eq.-1) xnr=cl/(-bl-al*dl2)
                !cc        xnro=xnr
                !cc       end if
                !cc       return
                !cc      end if
                !!!!!!!!!!!!!!

                !!new variant:
                izn=1
                dl2=-dsqrt(dll)/al
                xnr=-bl/al+dl2
                call dhdomega(pa,ptet,xnr,yn2)
                !cc        write(*,*)'#1: izn=',izn,' dl2=',dl2,' xnr=',xnr
                !cc        write(*,*)'znak=',znakstart,' -znak*dhdnr=',-znakstart*dhdnr
                if(-znakstart*dhdnr.gt.zero) then
                    izn=-1
                    dl2=dsqrt(dll)/al
                    xnr=cl/(-bl-al*dl2)
                    call dhdomega(pa,ptet,xnr,yn2)
                    !cc         write(*,*)'#2: izn=',izn,' dl2=',dl2,' xnr=',xnr
                    !cc         write(*,*)'znak=',znakstart,' -znak*dhdnr=',-znakstart*dhdnr
                    if(-znakstart*dhdnr.gt.zero) then
                        write(*,*)'Exception: both modes go outward !!'
                        stop
                    end if
                end if
                xnro=xnr
                !cc        pause
            end if
            return
        end if
        if(dll.lt.zero) goto(70,70,50) iroot
40      dl2=-dfloat(izn)*dsqrt(dll)/al
        if(izn.eq.1) xnr=-bl/al+dl2
        if(izn.eq.-1) xnr=cl/(-bl-al*dl2)
        xnro=xnr
        if(ivar.gt.1) then
            !cccccc  find Nr of reflected wave
            dnx=two*as*ynpopq+bs
            dhdnr=dnx*(two*g22*xnr-two*g12*yn2)/xj
            if(-znakstart*dhdnr.gt.zero) then
                izn=-izn
                goto 40
            end if
            return
        end if
50      if(iroot.eq.3) then
            !---------------------------
            !  find all roots
            !----------------------------
            if(dll.ge.zero) then
                xnr1=xnr
                xnr2=-bl/al-dl2
            else
                xnr1=1d+10
                xnr2=1d+10
            end if
            dll1=bl**2-al*cl1
            if(dll1.lt.zero) then
                xnr3=1d+10
                xnr4=1d+10
            else
                xnr3=-bl/al-izn*dsqrt(dll1)/al
                xnr4=-bl/al+izn*dsqrt(dll1)/al
            end if
        end if
        if(ider.eq.0) then
            prt=0d0
            prm=0d0
            return
        end if
        !--------------------------------------
        !   calculation of derivatives
        !--------------------------------------
        g11t=two*(dxdr*dxdrdt+dzdr*dzdrdt)
        g22t=two*(dxdt*dxdtdt+dzdt*dzdtdt)
        g33t=two*x0*(-pa*sitet-two*xgm*sitet*cotet)
        g12t=dxdrdt*dxdt+dxdr*dxdtdt+dzdrdt*dzdt+dzdr*dzdtdt
        xjt=g11t*g22+g22t*g11-two*g12*g12t
        btt=-b_tor*(r0/rm)/x0**2*x0t
        g2jqt=(g22t/xj-g22/xj**2*xjt)/(g2jq*two)
        bpt=xmy*(g2jqt*g3v-.5d0*g2jq*g3v/g33*g33t)
        bat=one/b*(bp*bpt+bt*btt)
        sit=bpt/b-bp/b**2*bat
        cot=btt/b-bt/b**2*bat
        u1t=c1/ww*bat
        e1t=-v/u**2*two*u1*u1t
        e2t=-v/u*u1t
        ynzt=yn2*sit*g2v1-yn2*si*g2v1**3/two*g22t + yn3*cot*g3v-yn3*co*g3v**3/two*g33t
        p1=two*ynz*ynzt
        p2=e1t
        p3=(e2*e2t)/e1-e2**2/(two*e1**2)*e1t

        s1=-p2/(two*e1**2)*e3*(ynzq-e1) + (e3+e1)/(two*e1)*(p1-p2)+p3
        s2=two*e3/e1*(ynzq-e1)*(p1-p2) - p2*e3/(e1**2)*(ynzq-e1)**2-two*e3*p3
        dnm=two*ynz*si*g2v1
        v1=(e3+e1)/(two*e1)*dnm
        v2=two*e3/e1*(ynzq-e1)*dnm
        !-----------------------------------
        !est !sav2009
        if(inew.gt.0) then
            gprt=-c0**2/ww**2*fnr/u*u1t*xsz
            if(inew.eq.1) then
                ynyt= - (yn2*cot*g2v1-yn2*co*g2v1**3/two*g22t - (yn3*sit*g3v-yn3*si*g3v**3/two*g33t))
                dnym= - co*g2v1
            else if(inew.eq.2) then
                ynyt= - g2jq*(yn2*cot*g2v1-yn2*co*g2v1**3/two*g22t - (yn3*sit*g3v-yn3*si*g3v**3/two*g33t))
                ynyt=ynyt - g2jqt*(yn2*g2v1*co-yn3*g3v*si)
                dnym= - g2jq*co*g2v1
            end if
            pnewt=(ynyt*gpr+yny*gprt)
            s1=s1+pnewt/(two*e1)-pnew/(two*e1**2)*e1t
            s2=s2+(pnewt*(ynzq-e3)+pnew*p1)/e1-pnew*(ynzq-e3)/e1**2*e1t
            v1=v1+dnym*gpr/(two*e1)
            v2=v2+gpr*(dnym*(ynzq-e3)+yny*dnm)/e1
        end if
        !---------------------------------------------
        vvt=-s1+(bs/as*s1-s2)/(two*dl1)
        vvm=-v1+(bs/as*v1-v2)/(two*dl1)
        s1=-yn2*(g12t/g22-g12/g22**2*g22t)
        s21=yn2**2*(g11t/g22-g11/g22**2*g22t)
        s22=yn3**2*( xjt/(g33*g22)-xj/(g33*g22)**2*(g33t*g22+g22t*g33) )
        sjg=(xjt*g22-xj*g22t)/g22**2
        s23=two*ynz*ynzt*xj/g22+sjg*ynzq
        s24=vvt*xj/g22+ynpopq*sjg
        s2=s21+s22-s23-s24
        prt=-s1+(two*(bl/al)*s1-s2)/(two*dl2)
        s1=-g12/g22
        s21=two*yn2*g11/g22
        s22=dnm*xj/g22
        s23=vvm*xj/g22
        s2=s21-s22-s23
        prm=-s1+(two*(bl/al)*s1-s2)/(two*dl2)
        if(ipow.gt.0) then
            !--------------------------------------
            !  calculation of decrements
            !--------------------------------------
            dnx=two*as*ynpopq+bs
            dhdnr=dnx*(two*g22*xnr-two*g12*yn2)/xj
            sl1=(ynzq-e1)*(ynzq+ynpopq-e1)-e2**2
            cf3=ynz
            cf4=xnr
            cf5=yn2
            vz=cltn/dabs(ynz)
            if(vz.gt.cltn) vz=cltn !sav2010
            vt=fvt(pa)
            jr=jfoundr
            icf1=iw
            icf2=izn
            call distr(vz,jr,ifound,fder)
            dfdv=fder
            vfound=vz
            cf2=ptet
            cf6=yny
            aimh=wpq/ww**2*pi*sl1*cltn**2/ynzq
            pdecv=dabs(aimh/dhdnr/xsz)
            !!        pdec1=-pdecv*dfdv
            pdec1=dabs(pdecv*dfdv)
            pnye=cnye*wpq**2/(pn*vt**3)
            pnyi=cnyi*pnye*zefff(pa)
            pdec2=dabs(pnyi/ww*(wpq/whe**2*ynpopq+wpq/ww**2*ynzq)*ynpopq/dhdnr/xsz)
            cf1=dsqrt(ynpopq)
            if(itend0.gt.0) then
                tmp=ft(pa)/0.16d-8
                fcoll=.5d-13*pn*zalfa**2*xlog/xmalfa/tmp**1.5d0
                !cc          ddens=dn1*pn
                !cc          tdens=dn2*pn
                !cc          tt=fti(pa)**0.33333d0    ! (ti, kev)^1/3
                !cc          source=4d-12*factor*ddens*tdens*dexp(-20d0/tt)/tt**2
                call source_new(pa,source)
                dek1=cnstal*pdecv*(1.d0-e3/ynpopq)**2/cf1
                dek2=source/(fcoll*pn)
                pdecal=dek1*dek2
                pdec3=zero
                if(itend0.gt.0) then
                    argum=clt/(cf1*valfa)
                    dek3=zatukh(argum,jr,vperp,kv)
                    pdec3=pdecal*dek3
                end if
            end if
        end if
        return
    !   conversion
    60    iconv=1
        if (ivar.ne.0) ivar=-1
        return
    !    reflection
    70    irefl=1
        if (ivar.gt.1.and.ivar.ne.10) then
            iw=-iw
            ivar=10
            goto 30
        end if
        if (ivar.eq.10) ivar=-1
        return
    end

    subroutine dhdomega(rho,theta,yn1,yn2)
        !use dispersion_module, only: yn3            
        implicit real*8 (a-h,o-z)
        !common /a0ef2/ ww
        !common /abefo/ yn3
        common/fj/dhdm,dhdnr,dhdtet,dhdr,ddn,dhdn3,dhdv2v,dhdu2u
        common/direct/znakstart
        parameter(zero=0.d0,h=1.d-6)
  
        call disp4(rho,theta,yn1,yn2)
  
  !!w*dH/dw=wdhdw:
        wdhdw=-(yn1*dhdnr+yn2*dhdm+yn3*dhdn3+dhdv2v+dhdu2u)
        znak=dsign(1.d0,wdhdw)
        znakstart=znak
  !c      write(*,*)'formula: znak=',znak
  !c      write(*,*)'wdhdw=',wdhdw,' H=',ham
  !c      write(*,*)'rho=',rho,' teta=',theta
  !c      write(*,*)'yn1=',yn1,' yn2=',yn2
  !c      write(*,*)'dhdnr=',dhdnr,' dhdm=',dhdm
  !c      write(*,*)'dhdr=',dhdr,' dhdtet=',dhdtet
  !c      write(*,*)'dhdn3=',dhdn3,' yn3=',yn3
  !c      write(*,*)'yn1*dhdnr=',yn1*dhdnr,' yn2*dhdm=',yn2*dhdm
  !c      write(*,*)'yn1*dhdnr+yn2*dhdm=',yn1*dhdnr+yn2*dhdm
  !cc      pause
  
        end    

    subroutine source_new(r,out)
        implicit real*8 (a-h,o-z)
        integer npta, klo, khi, ierr
        common /asou/ rsou(102),sou(102),npta
        call lock2(rsou,npta,r,klo,khi,ierr)
        if(ierr.ne.0) then
            write(*,*)'lock2 error in source_new'
            write(*,*)'ierr=',ierr,' rho=',r
            stop
        else
            call linf(rsou,sou,r,fout,klo,khi)
            out=dabs(fout)
        end if
    end        
end module dispersion_module