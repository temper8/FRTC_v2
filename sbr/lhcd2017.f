      subroutine lhcd2017(outpe)
cc******************************************************************
cc   outj(i)  = LH driven current density, MA/m^2
cc   outpe(i) =LH power density (Landau+coll.) deposited into electrons, MW/m^3
cc   outpec(i) = LH power density (collisions) deposited into electrons, MW/m^3
cc   outpef(i) = LH power dens. dep. into el. by fast wave, MW/m^3
cc   dndt(i)  = d^2Jr1/dt^2/E, MA/m^2/sec^2/(V/m), ~runaway d(el.density)/dt/E
cc   djdt(i)  = dJr2/dt, time drivative of runaway current Jr2, MA/m^2/sec
cc   outpa(i)  = LH power density deposited into alphas, MW/m^3
cc   outda(i)  = (Na/Ne) relative density of alpha-particles
cc******************************************************************
      use approximation
      use plasma
      use rt_parameters
      use spectrum1D      
      use maxwell  
      use spectrum_mod    
      use lhcd_module  
      implicit none
      integer i
      real*8 p_in,pe_p,pe_m,c_p,c_m
      real*8 vint
      include 'for/parameter.inc'
      include 'for/const.inc'
      include 'for/status.inc'
      real*8 outpe(NRD)
      real*8,dimension(:),allocatable:: outpep,outpem
      type(spectrum) spectr

cc*********************************************************************
cc    Co-ordinates used in ray-tracing:
cc         (x-x0)/rm=r*cos(teta)-delta-gamma*sin^2(teta) !sav2008 -gamma
cc         (z-z0)/rm=ell*r*sin(teta)
cc    Definitions:
cc    (x0,z0) - magnetic axis position, centimeters
cc    rm      - minor radius in mid-plane, cenrimeters
cc    r(rho_ASTRA),delta(r),gamma(r),ell(r) - dimensionless functions
cc    rho_ASTRA=sqrt(Phi_tor/GP/BTOR)
cc    Interval for r:  0.<= r <=1.
cc*********************************************************************
      print *, 'start start'
      tcur=time
      outpe=zero
      p_in=dble(QLH)    ! input LH power, MW

      if(p_in.eq.zero) then 
            dij(:,:,:)=zero
            return
      end if

      call read_parameters('lhcd/ray_tracing.dat')

      call init_plasma(NA1,ABC,BTOR,RTOR,UPDWN,GP2,
     & AMETR,RHO,SHIF,ELON,TRIA,MU,NE,TE,TI,ZEF,UPL)

      full_spectrum = read_spectrum('lhcd/spectrum.dat')
      full_spectrum%input_power = p_in
      pos_spectr = full_spectrum%get_positive_part()
      neg_spectr = full_spectrum%get_negative_part()

!!!!!!!!!!!!! starting ray-tracing !!!!!!!!!!!!!!!!!!!!!
      allocate(outpep(ngrid),outpem(ngrid))

!!positive spectrum:
      print *, 'positive spectrum'
      pe_p=zero
      outpep=zero
      if(pos_spectr%input_power > zero) then
            print *, 'spectrum_type',spectrum_type
            select case (spectrum_type)
            case (0)
                  spectr = make_spline_approximation(pos_spectr)
            case (1)
                  spectr = pos_spectr 
                  call spectr%calc_max_power
            case (2)
                  spectr = pos_spectr
                  call spectr%calc_max_power
            case (3)
                  print *, '2D spectrum'
                  stop                  
            end select
            call ourlhcd2017(spectr, outpep,pe_p)
      else
            dij(:,:,1)=zero
      end if      

      if(pe_p.ne.zero) then
            c_p=vint(outpep,roc)
            if(c_p.ne.zero) then
                  do i=1,ngrid
                        outpep(i)=pe_p*outpep(i)/c_p
                  end do
            end if
      end if

!!negative spectrum:
       print *, 'negative spectrum'
       pe_m=zero
       outpem=zero       
       if(neg_spectr%input_power > zero) then        
            select case (spectrum_type)
            case (0)
                  spectr = make_spline_approximation(neg_spectr)
            case (1)
                  spectr = neg_spectr 
                  call spectr%calc_max_power
            case (2)
                  spectr = neg_spectr
                  call spectr%calc_max_power
            case (3)
                  print *, '2D spectrum'
                  stop
            end select            
            call ourlhcd2017(spectr, outpem,pe_m)              
       else
            dij(:,:,2)=zero
       endif     
       
       if(pe_m.ne.zero) then
            c_m=vint(outpem,roc)
            if(c_m.ne.zero) then
                  do i=1,ngrid
                        outpem(i)=pe_m*outpem(i)/c_m
                  end do
            end if
      end if

      do i=1,ngrid
       outpe(i)=outpep(i)+outpem(i)
      end do

      deallocate(outpep,outpem)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dql1(pabs) !sav2008
      use rt_parameters
      use plasma, only : fvt, vperp
      use current
      !use spectrum1D, only: pabs
      use trajectory
      use dispersion_module
      use manager_mod, only: pow, inak, lenstor, lfree
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
      parameter(clt=3.d10,zero=0.d0)
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
c--------------------------------------
c   find power
c--------------------------------------
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
      call dfind(j,i,vz,powpr,pil,pic,pia,dfsr,dcv
     &                         ,refr,vlf,vrt,ifast)
c-----------------------------------
c      memorize trajectory
c----------------------------------
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine driver2(ystart,x1,x2,xsav,hmin,h1, pabs) !sav2008
      use constants
      use plasma
      use rt_parameters
      use manager_mod
      use trajectory
      use dispersion_module
      implicit real*8 (a-h,o-z)
      external extd2
      real*8 pabs
      !common /abc/ rzz,tetzz,xmzz,iznzz,iwzz,irszz
      !common /abcd/ irs
      !common /abcde/ izn!,iw
      !common /abcdg/ iabsorp
      !common /bcef/ ynz,ynpopq
      !common /bcg/ hrad
      !common /cefn/ iconv,irefl
      !common /ceg/ ipow,jfoundr
      common /cmn/ ind
      dimension ystart(2)
      parameter(nvar=2)
      dimension yscal(nvar),y(nvar),dydx(nvar),yold(nvar),dyold(nvar)
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
c-----------------------------
c            start moving
c-----------------------------
      do nstp=1,maxstep2
c---------------------------------------
c netpoint control
c---------------------------------------
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
c---------------------------------------------
c made step to nontransparent zone-return back
c----------------------------------------------
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
c--------------------------------------
c memorize step data
c--------------------------------------
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
c--------------------------------------
c choose step size
c--------------------------------------
       dst3=(x-xsav)*(x+h-xsav)
       if(dst3.lt.zero.and.irep.eq.0) h=xsav-x
       if(x.gt.rbord.and.h.gt.zero) then
        ind=2
        go to 20
       end if
10     dst1=(x-rbord)*(x+h-rbord)
       dst2=x*(x+h)
       if((dst1.lt.zero.and.irs.eq.-1).or.dst2.lt.zero) then
        h=h/2.d0
        if(dabs(h).lt.hmin1) then
         ind=4
         go to 20
        end if
        go to 10
       end if
c--------------------------------------
c find solution at x=x+hdid
c---------------------------------------
       ynz0=ynz
       call difeq(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,extd2)
20     continue
       if(ind.ne.0) then !exit
        xsav=xsav+hsav
        x2=x
        do i=1,nvar
         ystart(i)=y(i)
        end do
        ynz=ynz0
        return
       end if
c---------------------------------------
       if(dabs(hnext).lt.hmin) then
        if(ipri.gt.1) write(*,*) 'exit driver2: step is too small'
        go to 40
       end if
       h=hnext
      end do
c---------------------------------------
      if (ipri.gt.1) write (*,*) 'error in driver2: too many steps'
40    iabsorp=-1
      return
1001  format (10(e14.7,1x))
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine difeq(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs)
      use rt_parameters, only : hmin1
      implicit none
      external derivs
      integer nv,nmax,kmaxx,imax
      double precision eps,hdid,hnext,htry,x,dydx(nv),y(nv),yscal(nv)
     *,safe1,safe2,redmax,redmin,tiny,scalmx
     &,dysav(nv) !sav#
      parameter(nmax=50,kmaxx=8,imax=kmaxx+1,safe1=.25d0,safe2=.7d0
     *,redmax=1.d-5,redmin=.7d0,tiny=1.d-30,scalmx=.1d0)
cu    uses derivs,mmid,pzextr
      integer i,iq,k,kk,km,kmax,kopt,nseq(imax)
      double precision eps1,epsold,errmax,fact,h,red,scale,work,wrkmin
     *,xest,xnew,a(imax),alf(kmaxx,kmaxx),err(kmaxx)
     *,yerr(nmax),ysav(nmax),yseq(nmax)
      logical first,reduct
      save a,alf,epsold,first,kmax,kopt,nseq,xnew
      double precision dyd
      integer ii,ind
      common /cmn/ ind
      data first/.true./,epsold/-1.d0/
      data nseq /2,4,6,8,10,12,14,16,18/
      if(eps.ne.epsold)then
        hnext=-1.d29
        xnew=-1.d29
        eps1=safe1*eps
        a(1)=nseq(1)+1
        do k=1,kmaxx
          a(k+1)=a(k)+nseq(k+1)
        enddo
        do iq=2,kmaxx
          do k=1,iq-1
            alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.d0)*(2*k+
     *1)))
          enddo
        enddo
        epsold=eps
        do kopt=2,kmaxx-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1
        enddo
1      continue
       kmax=kopt
      endif
      h=htry
      do i=1,nv
        ysav(i)=y(i)
        dysav(i)=dydx(i)
      enddo
      if(h.ne.hnext.or.x.ne.xnew)then
        first=.true.
        kopt=kmax
      endif
      reduct=.false.
2     do k=1,kmax
        xnew=x+h
        if(xnew.eq.x) then
            write(*,*) 'step size underflow in difeq'
            pause
        end if
        call mmid(ysav,dydx,nv,x,h,nseq(k),yseq,derivs)
 !sav#
      if(ind.eq.1) then
        h=h/2d0
        if (dabs(h).lt.hmin1) then
          do ii=1,nv
            y(ii)=ysav(ii)
          end do
          hnext=h
          return
        end if
        do ii=1,nv
          dyd=dabs(dysav(ii))
          yscal(ii)=dabs(ysav(ii))+dabs(h*dyd)+1.d-30/(1d0+dyd)
          y(ii)=ysav(ii)
        end do
        goto 2
      end if
!sav#
        xest=(h/nseq(k))**2
!var        call pzextr(k,xest,yseq,y,yerr,nv)  !polynomial extrapolation
        call rzextr(k,xest,yseq,y,yerr,nv) !rational extrapolation
        if(k.ne.1)then
          errmax=tiny
          do i=1,nv
            errmax=max(errmax,abs(yerr(i)/yscal(i)))
          enddo
          errmax=errmax/eps
          km=k-1
          err(km)=(errmax/safe1)**(1.d0/(2*km+1))
        endif
        if(k.ne.1.and.(k.ge.kopt-1.or.first))then
          if(errmax.lt.1.d0)goto 4
          if(k.eq.kmax.or.k.eq.kopt+1)then
            red=safe2/err(km)
            goto 3
          else if(k.eq.kopt)then
            if(alf(kopt-1,kopt).lt.err(km))then
              red=1.d0/err(km)
              goto 3
            endif
          else if(kopt.eq.kmax)then
            if(alf(km,kmax-1).lt.err(km))then
              red=alf(km,kmax-1)*safe2/err(km)
              goto 3
            endif
          else if(alf(km,kopt).lt.err(km))then
            red=alf(km,kopt-1)/err(km)
            goto 3
          endif
        endif
      enddo
3     red=min(red,redmin)
      red=max(red,redmax)
      h=h*red
      reduct=.true.
      goto 2
4     x=xnew
      hdid=h
      first=.false.
      wrkmin=1.d35
      do kk=1,km
        fact=max(err(kk),scalmx)
        work=fact*a(kk+1)
        if(work.lt.wrkmin)then
          scale=fact
          wrkmin=work
          kopt=kk+1
        endif
      enddo
      hnext=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
        fact=max(scale/alf(kopt-1,kopt),scalmx)
        if(a(kopt+1)*fact.le.wrkmin)then
          hnext=h/fact
          kopt=kopt+1
        endif
      endif
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs)
      use dispersion_module, only: iconv,irefl
      implicit none
      external derivs
      integer nstep,nvar,nmax
      double precision htot,xs,dydx(nvar),y(nvar),yout(nvar)
      parameter (nmax=50)
      integer i,n
      double precision h,h2,swap,x,ym(nmax),yn(nmax)
      double precision yz1,yz2
      !integer iconv,irefl
      !common /cefn/ iconv,irefl
      integer ind
      common /cmn/ ind
      h=htot/nstep
      yz1=y(1) !sav#
      yz2=y(2) !sav#
      do 11 i=1,nvar
        ym(i)=y(i)
        yn(i)=y(i)+h*dydx(i)
11    continue
      x=xs+h
      call derivs(x,yn,yout)
      if (iconv+irefl.ne.0) goto 10 !sav#
      h2=2.d0*h
      do 13 n=2,nstep
        do 12 i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
12      continue
        x=x+h
        call derivs(x,yn,yout)
        if (iconv+irefl.ne.0) goto 10 !sav#
13    continue
      do 14 i=1,nvar
        yout(i)=0.5d0*(ym(i)+yn(i)+h*yout(i))
14    continue
      ind=0 !sav#
      return
10    ind=1 !sav#
      yout(1)=yz1 !sav#
      yout(2)=yz2 !sav#
      return !sav#
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rzextr(iest,xest,yest,yz,dy,nv)
      integer iest,nv,imax,nmax
      double precision xest,dy(nv),yest(nv),yz(nv)
      parameter (imax=13,nmax=50)
      integer j,k
      double precision b,b1,c,ddy,v,yy,d(nmax,imax),fx(imax),x(imax)
      save d,x
      x(iest)=xest
      if(iest.eq.1) then
        do 11 j=1,nv
          yz(j)=yest(j)
          d(j,1)=yest(j)
          dy(j)=yest(j)
11      continue
      else
        do 12 k=1,iest-1
          fx(k+1)=x(iest-k)/xest
12      continue
        do 14 j=1,nv
          yy=yest(j)
          v=d(j,1)
          c=yy
          d(j,1)=yy
          do 13 k=2,iest
            b1=fx(k)*v
            b=b1-c
            if(b.ne.0.d0) then
              b=(c-v)/b
              ddy=c*b
              c=b1*b
            else
              ddy=v
            endif
            if (k.ne.iest) v=d(j,k)
            d(j,k)=ddy
            yy=yy+ddy
13        continue
          dy(j)=ddy
          yz(j)=yy
14      continue
      endif
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine pzextr(iest,xest,yest,yz,dy,nv)
      integer iest,nv,imax,nmax
      double precision xest,dy(nv),yest(nv),yz(nv)
      parameter (imax=13,nmax=50)
      integer j,k1
      double precision delta,f1,f2,q,d(nmax),qcol(nmax,imax),x(imax)
      save qcol,x
      x(iest)=xest
      do 11 j=1,nv
        dy(j)=yest(j)
        yz(j)=yest(j)
11    continue
      if(iest.eq.1) then
        do 12 j=1,nv
          qcol(j,1)=yest(j)
12      continue
      else
        do 13 j=1,nv
          d(j)=yest(j)
13      continue
        do 15 k1=1,iest-1
          delta=1.d0/(x(iest-k1)-xest)
          f1=xest*delta
          f2=x(iest-k1)*delta
          do 14 j=1,nv
            q=qcol(j,k1)
            qcol(j,k1)=dy(j)
            delta=d(j)-q
            dy(j)=f1*delta
            d(j)=f2*delta
            yz(j)=yz(j)+dy(j)
14        continue
15      continue
        do 16 j=1,nv
          qcol(j,iest)=dy(j)
16      continue
      endif
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine distr(vz,j,ifound,fder)
      use iterator_mod
      use lock_module      
      implicit none
      real*8 vz,fder
      integer j,ifound,i,klo,khi,ierr,nvp
      real*8,dimension(:),allocatable:: vzj,dfdvj
      real*8 vlf,vrt,dflf,dfrt,dfout
      !common/gridv/vgrid(101,100),dfundv(101,100),nvpt
      common /a0ghp/ vlf,vrt,dflf,dfrt
      nvp=nvpt
      allocate(vzj(nvp),dfdvj(nvp))
      do i=1,nvp
       vzj(i)=vgrid(i,j)
       dfdvj(i)=dfundv(i,j)
      end do
      call lock2(vzj,nvp,vz,klo,khi,ierr)
      if(ierr.eq.0) then !vgrid(1,j) <= vz <= vgrid(nvpt,j)
       call linf(vzj,dfdvj,vz,dfout,klo,khi)
       ifound=klo
       vlf=vzj(klo)
       vrt=vzj(khi)
       fder=dfout
       dflf=dfdvj(klo)
       dfrt=dfdvj(khi)
      else if(ierr.eq.1) then !vz < vgrid(1,j)
       write(*,*)'exception: ierr=1 in distr()'
       pause'next key = stop'
       stop
      else if(ierr.eq.2) then !vz > vgrid(nvpt,j)
       write(*,*)'exception: ierr=2 in distr()'
       pause'next key = stop'
       stop
      else if(ierr.eq.3) then
       write(*,*)'exception in distr, klo=khi=',klo,' j=',j,' nvp=',nvp
       write(*,*)'vz=',vz,' v1=',vzj(1),' v2=',vzj(nvp)
       pause'next key = stop'
       stop
      end if
      deallocate(vzj,dfdvj)
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine interpol(fj)
      implicit none
      integer p,k,i0,i
      real*8 p2,zero
      double precision fj(1002)
      parameter(zero=0.d0,i0=1002)
      p=zero
      do i=1,i0
            if (fj(i).gt.zero) then
                  if (p.gt.zero) then
                        p2=(fj(i)-fj(p))/(i-p)
                              do k=p+1,i-1
                                    fj(k)=fj(p)+p2*(k-p)
                              end do
                        end if
                  p=i
                  end if
            end do
      end
