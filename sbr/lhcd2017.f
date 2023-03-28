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
      subroutine ourlhcd2017(spectr, outpe,pe_out)      
      use constants
      use approximation
      use spline_module
      use chebyshev
      use plasma
      use rt_parameters
      !use spectrum1D
      use maxwell      
      use trajectory, only: view, nrefj, init_trajectory
      use spectrum_mod
      use manager_mod
      use dispersion_module
      use current
      use iteration_result_mod
      use iterator_mod
      use lock_module
      implicit real*8 (a-h,o-z)
      type(spectrum) spectr
      real*8 outpe,pe_out 
      dimension outpe(*)
      dimension galfa(50,100),vpmin(100),vcva(100)
     &,pd2(100),pd2a(100),pd2b(100),pdprev1(100),pdprev2(100)
     &,source(100),sour(100)
     &,rxx(102),pwe(102),wrk(102)
      !dimension vmid(100),vz1(100),vz2(100),ibeg(100),iend(100)
      !common /a0a4/ plost,pnab
      !common /bcef/ ynz,ynpopq
      common /a0ghp/ vlf,vrt,dflf,dfrt
      common/plosh/ zv1(100,2),zv2(100,2)!,sk(100)
      !common /asou/ rsou(102),sou(102),npta
      !common/gridv/vgrid(101,100),dfundv(101,100),nvpt
      !common /vvv2/ psum4
      !common /arr/ dgdu(50,100),kzero(100)
      !common /ag/ inak,lenstor,lfree
      common /maxrho/ rmx_n,rmx_t,rmx_z,rmx_ti
      
      type(IterationResult) :: iteration_result
      real*8 kofpar,timecof
      !real*8,dimension(:),allocatable:: vvj,vdfj

      !double precision vrj(101),dj(101),djnew(1001)
      !double precision dj2(101),d2j(101)

      integer iptnew
      real*8 dijk, vrjnew, plaun
      common/t01/dijk(101,100,2), vrjnew(101,100,2), iptnew
      integer ispectr
      plaun = spectr%input_power

      ispectr = spectr%direction
      lfree=1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      hr = 1.d0/dble(nr+1)
      iw0=iw
      nrr=nr+2
      rxx(1)=zero
      rxx(nrr)=one
      do j=1,nr
            rxx(j+1)=hr*dble(j)
      end do
      
      call find_volums_and_surfaces

      ppv1=zero
      ppv2=zero
      pnab=zero
      plost=zero
      psum4=zero
      anb=zero
      fuspow=zero
      o_da=zero
c-------------------------------------------
c find velocity limits and initial dfdv
c--------------------------------------------
      ipt1=kpt1+1
      ipt2=ni1+ni2
      ipt=ipt1+ni1+ni2+kpt3
      if(ipt.gt.101) then
            write(*,*)'ipt >101'
            pause'stop program'
            stop
      end if
      nvpt=ipt

      do j=1,nr                  ! begin 'rho' cycle
            r=hr*dble(j)
!!!!sav2008       pn=fn(r)
!!       pn=fn1(r,fnr)
!!       pn=fn2(r,fnr,fnrr) !sav2008
            if(inew.eq.0) then !vardens
                  pn=fn1(r,fnr)
            else
                  pn=fn2(r,fnr,fnrr)
            end if
            dens(j)=pn
            vt=fvt(r)
            vto=vt/vt0
            wpq=c0**2*pn
            whe=dabs(b_tor)*c1
            v=wpq/ww**2
            u1=whe/ww
            u=u1**2
            e1=1d0-v*(1d0/xmi-1d0/u)
            e2=v/u1
            e3=v
            tmp=ft(r)/0.16d-8 !Te, keV
            cn1=dsqrt(50d0/tmp)  !sav2008
            if(itend0.gt.0) then
                  eta(j)=1d0-v
                  vcva(j)=cnstvc*vt*dsqrt(2d0)/valfa
                  vpmin(j)=2.0d0*dsqrt(tmp/(-eta(j)))
222               continue
                  dvperp=(vpmax-vpmin(j))/dble(kv-1)
                  if(dvperp.le.zero) then
                        vpmax=1.3d0*vpmax
                        go to 222
                  end if
                  do k=1,kv
                        vperp(k,j)=vpmin(j)+dble(k-1)*dvperp
                  end do
               fcoll(j)=.5d-13*dens(j)*zalfa**2*xlog/xmalfa/tmp**1.5d0
                  ddens=dn1*dens(j)
                  tdens=dn2*dens(j)
                  tt=fti(r)**one_third    ! (ti, keV)^1/3
               source(j)=4d-12*factor*ddens*tdens*dexp(-20d0/tt)/tt**2
                  anb=anb+source(j)*vk(j)
            end if
            cn2=dsqrt(dabs(e1))+e2/dsqrt(e3) !sav2008
            !vz1(j)=cleft*cltn/cn1  !Vpar/Vt0
            !vz2(j)=cright*cltn/cn2  !Vpar/Vt0
            !if(vz2(j).gt.0.9d0*cltn) vz2(j)=0.9d0*cltn
            !v1=vz1(j)/vto !Vpar/Vt(rho)
            !v2=vz2(j)/vto !Vpar/Vt(rho)
            vmax=cltn/vto
            v1=4.d0  !Vpar/Vt(rho)
            v2=10.d0 !cright*cltn/cn2 !10.d0 !Vpar/Vt(rho)
            if(v2.ge.vmax) v2=0.5d0*vmax
            if(v1.ge.v2) v1=v2-2.d0
            call gridvel(v1,v2,vmax,0.5d0,ni1,ni2,ipt1,kpt3,vrj)
            vz1(j)=v1*vto !Vpar/Vt0
            vz2(j)=v2*vto !Vpar/Vt0
            if(vz2(j).gt.0.9d0*cltn) vz2(j)=0.9d0*cltn
            do i=1,ipt
                  vgrid(i,j)=vrj(i)*vto
            end do
      end do                     ! end 'rho' cycle 

!!!!!!!!!read data !!!!!!!!!!!!       
      allocate(vvj(i0),vdfj(i0))
      k=(3-ispectr)/2
      do j=1,nr
            r=hr*dble(j)
            vt=fvt(r)
            vto=vt/vt0
            do i=1,i0
                  vvj(i)=vij(i,j)
                  vdfj(i)=dfij(i,j,k) !=dfundv(i,j)*vto**2
            end do
            do i=1,ipt
                  vrj(i)=vgrid(i,j)/vto   !Vpar/Vt
                  call lock(vvj,i0,vrj(i),klo,khi,ierr)
            if(ierr.eq.1) then
                  write(*,*)'lock error in read distribution function'
                  write(*,*)'j=',j,'i0=',i0
                  write(*,*)'vvj(1)=',vvj(1),' vvj(i0)=',vvj(i0)
                  write(*,*)'i=',i,' vrj(i)=',vrj(i),' vmax=',cltn/vto
                  write(*,*)
                  pause'next key = stop'
                  stop
            end if
            call linf(vvj,vdfj,vrj(i),dfout,klo,khi)
            dfundv(i,j)=dfout/vto**2
            if(dfundv(i,j).gt.zero) dfundv(i,j)=zero
            end do
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(itend0.gt.0) then  ! begin alpha-source renormalisation
            fuspow=anb*talfa*1.6022d-19
            anb0=anb
            anb=zero
            do j=1,nr
                  r=hr*dble(j)
                  if(r.le.dra) then
                        tt=fti(zero)**one_third
                  else
                  tt=fti(r-dra)**one_third    ! (shifted ti, kev)^1/3
                  end if
                  ddens=dn1*dens(j)
                  tdens=dn2*dens(j)
                  sour(j)=4d-12*factor*ddens*tdens*dexp(-20d0/tt)/tt**2
                  anb=anb+sour(j)*vk(j)
            end do
            aratio=anb0/anb
            rsou(1)=zero
            sou(1)=aratio*sour(1)
            do j=1,nr
                  r=hr*dble(j)
                  rsou(j+1)=r
                  sou(j+1)=aratio*sour(j)
                  if(j.eq.nr) sssour=source(j)
                  source(j)=sou(j+1)
            end do
            npta=nr+2
            rsou(npta)=1.d0
            sou(npta)=aratio*sour(nr)
      end if
c------------------------------------
c set initial values of arrays
c------------------------------------
      dland=zero
      dcoll=zero
      perpn=zero
      dalf=zero
      vel=zero
      jrad=0
      iww=0
      tatai=zero
      xnpar=zero
      izz=zero
!
      pdl=zero
      pdc=zero
      pda=zero
      pdfast=zero
      pdprev1=zero
      pdprev2=zero
      tok=zero
      cur=zero
      pd2=zero
      pd2a=zero
      pd2b=zero
      dql=zero
      dq1=zero
      dq2=zero
      dncount=zero
      vzmin=cltn
      vzmax=-cltn
      kzero=kv
      call init_trajectory

      if(itend0.gt.0) then
            do j=1,nr           ! begin 'rho' cycle
                  do i=1,50
                        dqi0(i,j)=zero
                  end do
                  call alphas(dqi0,vperp,j,kv,galfa)
            end do              ! end 'rho' cycle
      end if
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!sign of driven current in right coordinate system {dro,dteta,dfi}:
!!!!!curdir=+1.0 for current drive in positive direction "dfi"
!!!!!curdir=-1.0 for current drive in negative direction "dfi"
!!!!!spectrum Nz>0 is along dfi>0 and Nz<0 is along dfi<0
!!!!!it is also OK if Npar is used instead of Nz, but for Btor>0, that is along dfi>0
!!      curdir=-dble(ispectr)
!!!!!!!!!!!!!!! begin iterations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      q_rest=plaun
      iterat=0

80    continue
      call manager(iterat, iw0, ntet, spectr)

      call find_achieved_radial_points(nvpt)

      do j=1,nr
            pdl(j)=pdl(j)*xwtt
            pdc(j)=pdc(j)*xwtt
            pda(j)=pda(j)*xwtt
            pdfast(j)=pdfast(j)*xwtt
            pwe(j+1)=(pdl(j)+pdc(j))/vk(j)
      end do
      pwe(1)=pwe(2)
      pwe(nr+2)=zero

!!   find nevyazka
!!----------------------------
      psum1=zero
      psum2=zero
      pchg=zero
      pchg1=zero
      pchg2=zero
      do j=1,nr
            dpw1=pdl(j)+pdc(j)
            dpw2=pda(j)
            psum1=psum1+dpw1**2
            psum2=psum2+dpw2**2
            pchg1=pchg1+(dpw1-pdprev1(j))**2
            pchg2=pchg2+(dpw2-pdprev2(j))**2
            pdprev1(j)=dpw1
            pdprev2(j)=dpw2
      end do
      if(psum1.ne.zero) pchg=pchg1/psum1 !sav2008
      if(psum2.ne.zero) pchg=pchg+pchg2/psum2
c----------------------------------------
c     calculate total current and power
c----------------------------------------
      cppl=zero
      cppc=zero
      cppa=zero
      cppf=zero
      do j=1,nr
            cppl=cppl+pdl(j)
            cppc=cppc+pdc(j)
            cppa=cppa+pda(j)
            cppf=cppf+pdfast(j)
      end do
      ol=cppl*1d-6
      oc=cppc*1d-6
      oa=cppa*1d-6
      of=cppf*1d-6
!!!!!!!!! prepare to the next iteration !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      iterat=iterat+1
      q_abs=ol+oc+oa
      q_rest=pnab/xsgs
      q_cond=zero
      if(q_abs.ne.zero) q_cond=0.5d0*q_rest/q_abs
!!!      if(q_cond.le.pabs0.and.pchg.lt.pgiter)

      call integral(1,nspl,rh,con,avedens) 

      iteration_result = IterationResult(number = iterat, 
     & spectr_direction = ispectr, P_launched = plaun,
     & P_landau = ol, P_coll = oc, P_alph = oa,
     & alphas_power = fuspow, P_fast = of,
     & P_lost = plost/xsgs, P_not_accounted = pnab/xsgs,
     & P_landau_strong_absorption = ppv1/xsgs,
     & P_landau_weak_absorption = ppv2/xsgs,
     & P_turns = psum4/xsgs, efficiency = oi/plaun,
     & avedens = avedens*1.d19, r0 = r0*1.d-2,
     & eta_eff = 1.d17*avedens*r0*oi/plaun,
     & residual = pchg)

      if(iterat.gt.5.and.q_cond.le.pabs0.and.pchg.lt.pgiter) goto 110

      if(ipri.gt.1) then
            call iteration_result%print
            call iteration_result%save(tcur)
            !pause
      end if

      if(iterat.le.niterat) then
            call recalculate_f_for_a_new_mesh(ispectr)
            call init_iteration
            goto 80
      end if
c------------------------------------------
c save results
c------------------------------------------
110   continue

      if(ipri.gt.0) then
       write (*,*)
       write (*,*) 'RAY-TRACING RESULTS:'
       call iteration_result%print
       call iteration_result%save(tcur)
       write (*,*) '-------------------------------------------'
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      k=(3-ispectr)/2
      do j=1,nr
       r=hr*dble(j)
       vt=fvt(r)
       vto=vt/vt0
       vmax=cltn/vto
       zff=(5d0+zefff(r))/5d0
       cnyfoc=zff*c0**4*cnye
       if(inew.eq.0) then !vardens
        pn=fn1(r,fnr)
       else
        pn=fn2(r,fnr,fnrr)
       end if
       dconst=vt0/(1.d-10*cnyfoc*pme*pn**2) !divided by 10^-10 here 
!!!!!!!!                         and multiplied by 10^-10 in dfind()
!!!old       dconst=vt0/(cnyfoc*pme*pn**2)
!!!        dj(i)=dql(i,j)*dconst*vto !D_normir
       do i=1,ipt
        vrj(i)=vgrid(i,j)/vto      !Vpar/Vt
        dj(i)=dql(i,j)*dconst*vto  !D_normir
        vrjnew(i,j,k)=vrj(i)
        dijk(i,j,k)=dj(i)
       end do
       do i=1,i0
        if(vij(i,j).ge.vmax) then
         ddout=zero
        else
         call lock(vrj,ipt,vij(i,j),klo,khi,ierr)
         if(ierr.eq.1) then
          write(*,*)'lock error in output dql'
          write(*,*)'j=',j,'ipt=',ipt
          write(*,*)'vrj(1)=',vrj(1),' vrj(ipt)=',vrj(ipt)
          write(*,*)'i=',i,' v=',vij(i,j),' vmax=',vmax
          write(*,*)
          pause'next key = stop'
          stop
         end if
!!!         call linf(vrj,dj,vij(i,j),ddout,klo,khi)
!!         if(ddout.le.1.d0) ddout=zero
         ddout=dj(klo)
        end if
        dij(i,j,k)=ddout
       end do
       zv1(j,k)=vrj(ipt1)
       zv2(j,k)=vrj(ni1+ni2+ipt1)
      end do

      if (ispectr == -1) call view(tcur,1,spectr%size,ntet)  !writing trajectories into a file

      if(ismthout.ne.0) then
       do i=1,nrr
        wrk(i)=pwe(i)
       end do
       call fsmoth4(rxx,wrk,nrr,pwe)
      end if
!
      rh(1)=rh1
      if(rh(nspl).gt.1.d0) rh(nspl)=1.d0
      do j=1,nspl
       call lock2(rxx,nrr,rh(j),klo,khi,ierr)
       if(ierr.ne.0) then
        write(*,*)'lock2 error in profiles for ASTRA'
        write(*,*)'ierr=',ierr,' j=',j,' rh(j)=',rh(j)
        write(*,*)'rxx(1)=',rxx(1),' rxx(nrr)=',rxx(nrr)
        pause
       end if
       call linf(rxx,pwe,rh(j),fout,klo,khi)
       outpe(j)=fout
      end do
      pe_out=ol+oc
      rh(1)=zero
!
      deallocate(vvj,vdfj)
      end

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
      subroutine dqliter(dltpow,ib,ie,h,powexit,iout) !sav2008
      use rt_parameters
      use trajectory
      use dispersion_module
      use current
      use plasma, only: vperp
      use iterator_mod, only: psum4
      implicit real*8 (a-h,o-z)
      dimension an1(length),an2(length)
      common /xn1xn2/ an1,an2
      common /a0ghp/ vlf,vrt,dflf,dfrt
      !common /vvv2/ psum4
      parameter(clt=3.d10,zero=0.d0)
      pow=powexit
      pdec1=zero
      pdec1z=zero
      pdec3=zero
      pdec3z=zero
      pdecv=zero
      pintld=zero
      pintal=zero
10    continue
      iout=0
      do i=ib,ie
c-----------------------------------
c restore memorized decrements and
c integrate power equation
c------------------------------------
       v=vel(i)
       jr=jrad(i)
       refr=perpn(i)
       ifast=iww(i)
       dek3=zero
       if(itend0.gt.0) then
        argum=clt/(refr*valfa)
        dek3=zatukh(argum,abs(jr),vperp,kv)
       end if
!!!old variant
!!!       call raspr(v,abs(jr),iv,df)
!!!       if(iv.eq.0) iv=1
!!!!!!!!!!!!!!!!!!!!!!!!!!
       call distr(v,abs(jr),iv,df)
!!       dfsr=v*df*(vrt-vlf)
!!       vsr=v*(vrt-vlf)
       dfsr=(vlf*dflf+vrt*dfrt)/2d0*(vrt-vlf) !sav2008
       vsr=(vrt+vlf)*(vrt-vlf)/2d0 !sav2008
       if(jr.lt.0) then !case of turn
        jr=-jr
!variant        pintld=-dland(i)*df
!!        pintld=-dland(i)*(dflf+dfrt)/2d0
        pintld=dabs(dland(i)*(dflf+dfrt)/2d0)
        pdec2=dexp(-2d0*dcoll(i))
        pintal=dabs(dalf(i)*dek3)
        pcurr=pdec2*dexp(-2d0*pintld-2d0*pintal)
        psum4=psum4+pow*(1d0-pcurr)
        dcv=dland(i)/vsr
       else
        pdec2=dcoll(i)
        pdecv=dland(i)
!!        pdec1=-pdecv*df
        pdec1=dabs(pdecv*df)
        pdec3=dabs(dalf(i)*dek3)
        pintld=(pdec1+pdec1z)/2d0*h
        pintal=(pdec3+pdec3z)/2d0*h
        pdec1z=pdec1
        pdec3z=pdec3
        dcv=pdecv*h/vsr
       end if
       powpr=pow
       if(dltpow.ne.zero) then
        powd=pow*dexp(-2d0*pintld)
        powcol=powd*pdec2
        powal=powcol*dexp(-2d0*pintal)
        pow=powal
       end if
       pil=pintld
       pic=.5d0*dabs(dlog(pdec2))
       pia=pintal
       call dfind(jr,iv,v,powpr,pil,pic,pia,dfsr,dcv
     &                           ,refr,vlf,vrt,ifast)
       if(pow.lt.dltpow) then
        powexit=pow
        return
       end if
      end do
      jchek=jrad(ie+1)
c-------------------------------------------
c  check whether trajectory has continuation
c---------------------------------------------
      if(jchek.eq.0) then
       iout=1
       powexit=pow
       return
      else
       ib=idnint(dland(ie+1))
       ie=idnint(dcoll(ie+1))
       goto 10
      end if
      end
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
      subroutine driver4(ystart,x1,x2,rexi,hmin,derivs)
      use constants
      use Runge_Kutta_module
      use plasma
      use rt_parameters
      use trajectory, only: irs, iabsorp
      use dispersion_module
      implicit real*8 (a-h,o-z)
      external derivs
      !common /abcd/ irs
      !common /abcde/ izn!,iw
      !common /abcdg/ iabsorp
      !common /bdeo/ ivar
      !common /bcef/ ynz,ynpopq
      !common /df/ pdec14,pdec24,pdec34,idec
      common /dg/ pintld4,pintcl4,pintal4
      parameter(hbeg=1.d-4,iturns=1,maxat=3,nvar=4) !sav2008
      dimension ystart(4),yscal(nvar),y(nvar),dydx(nvar),yold(nvar)
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
10    continue
c--------------------------------------
c start integration
c--------------------------------------
      do nstp=1,maxstep4
        idec=iturns
        call derivs(x,y,dydx)
        idec=0
        pintld4=pintld4+dabs((pdec14+pdec14zz)/2d0*hdid)
        pintcl4=pintcl4+dabs((pdec24+pdec24zz)/2d0*hdid)
        pintal4=pintal4+dabs((pdec34+pdec34zz)/2d0*hdid)
        pdec14zz=pdec14
        pdec24zz=pdec24
        pdec34zz=pdec34
       if(nstp.eq.1) then
        h=hbeg
!!var        if(dabs(dydx(3)).ne.zero) h=dabs(hmin1/dydx(3))/hdrob1
        if(dabs(dydx(3)).ne.zero) h=0.5d0*dabs(rrange/dydx(3))/hdrob1
       end if
20     continue
       if(y(3).ge.rbord1.and.dydx(3).gt.zero) then
c--------------------------------------
c forced reflection from periphery
c--------------------------------------
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
1      continue
c---------------------------------------
c remember old values
c---------------------------------------
       xold=x
       do i=1,nvar
        dyd=dabs(dydx(i))
        yscal(i)=dabs(y(i))+dabs(h*dyd)+1.d-30/(1.d0+dyd)+1.d-30
        yold(i)=y(i)
       end do
       if (y(3)*irs.lt.rmm*irs) rmm=y(3)
30     continue
!!var       call rkqc(y,dydx,nvar,x,h,eps1,yscal,hdid,hnext,derivs)
       call rkqs(y,dydx,nvar,x,h,eps1,yscal,hdid,hnext,derivs)
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
      subroutine alphas(d,u,j,kmax,g)
      use dispersion_module, only: dgdu, kzero
      implicit real*8 (a-h,o-z)
      dimension d(50,100),u(50,100),g(50,100)
      !common /arr/ dgdu(50,100),kzero(100)
      parameter(zero=0.d0, one=1.d0, tiny=1.d-30)
        km=kzero(j)
        um=u(km,j)
        if(um.ge.one) then
          do k=1,kmax
            if(u(k,j).lt.one) then
             uk=u(k,j)
             uk2=uk**2
             w=dsqrt(one-uk2)
             g(k,j)=w/uk2
             dgdu(k,j)=-one/(w*uk)-2.d0*w/(uk*uk2)
            else
             g(k,j)=zero
             dgdu(k,j)=zero
            end if
          end do
          return
        end if

          do k=1,km
           uk=u(k,j)
           uk2=uk**2
           w=dsqrt(one-uk2)
           g(k,j)=w/uk2
           dgdu(k,j)=-one/(w*uk)-2.d0*w/(uk*uk2)
          end do

       do k=km+1,kmax
         du=u(k,j)-u(k-1,j)
          if(u(k,j).lt.one) then
            beta=u(k,j)*dsqrt(one-u(k,j)**2)
          else
            beta=zero
          end if
         alfa=u(k,j)**3
         g(k,j)=(d(k,j)*g(k-1,j)+beta*du)/(d(k,j)+alfa*du)
           if(d(k,j).ne.zero) then
            dgdu(k,j)=(beta-alfa*g(k,j))/d(k,j)
           else
            dgdu(k,j)=(g(k,j)-g(k-1,j))/du
           end if
       end do
               do k=1,kmax
                if(g(k,j).lt.tiny) g(k,j)=zero
                if(dabs(dgdu(k,j)).lt.tiny) dgdu(k,j)=zero
               end do
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine diff(x,y,n,dy)
      implicit real*8 (a-h,o-z)
      dimension y(*),x(*),dy(*)
        dy(1)=(y(2)-y(1))/(x(2)-x(1))
        do k=2,n-1
         dy(k)=(y(k+1)-y(k-1))/(x(k+1)-x(k-1))
        end do
        dy(n)=(y(n)-y(n-1))/(x(n)-x(n-1))
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine integral(ibeg,iend,x,y,fout)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      fout=0.d0
      if(ibeg.eq.iend) return
      znak=1.d0
      n1=ibeg
      n2=iend
      if(n2.lt.n1) then
       znak=-1.d0
       ie=n1
       n1=n2
       n2=ie
      end if
      sum=0.d0
      do i=n1+1,n2
       dx=x(i)-x(i-1)
       dsum=y(i)+y(i-1)
       sum=sum+.5d0*dsum*dx
      end do
      fout=znak*sum
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine gridvel(v1,v2,vmax,cdel,ni1,ni2,ipt1,kpt3,vrj)
      implicit none
      integer ni1,ni2,ipt1,kpt1,kpt2,kpt3,k
      double precision vrj(*),v1,v2,v12,vmax,cdel
      kpt1=ipt1-1
      kpt2=ni1+ni2+1
      do k=1,kpt1  !0<=v<v1
       vrj(k)=dble(k-1)*v1/dble(kpt1)
      end do
      v12=v1+(v2-v1)*cdel
      do k=1,ni1+1 !v1<=v<=v12
       vrj(k+kpt1)=v1+dble(k-1)*(v12-v1)/dble(ni1)
      end do
      do k=2,ni2+1 !!v12<v<=v2
       vrj(k+kpt1+ni1)=v12+dble(k-1)*(v2-v12)/dble(ni2)
      end do     
      do k=1,kpt3  !v2<v<=vmax
       vrj(k+kpt1+kpt2)=v2+dble(k)*(vmax-v2)/dble(kpt3)
      end do
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fsmoth4(x,y,n,ys)
      use approximation
      implicit real*8 (a-h,o-z)
      !external polin2
      parameter(zero=0.d0,np=10,imax=601)
      parameter(m0=1,ndp=1)
      ! m0,ndp - parameters of smoothing procedure
      dimension y(n),x(n),ys(n)
      dimension yy(imax),xx(imax)
      dimension coeffs(np),cffs(np)
      dimension dys(imax)
      if(n.gt.imax) stop 'small imax in subroutine fsmoth4()'
        call diff(x,y,n,dys)
         do k=1,n
          ys(k)=y(k)
         end do
       m=m0
       m2=m+2
       id=m+ndp
       nmax=n-id
       xs=x(1)
       do j=1,nmax
         do i=1,id
          xx(i)=x(j+i)-xs
          yy(i)=y(j+i)-ys(j)-dys(j)*xx(i)
         end do
        call approx(xx,yy,id,polin2,m,coeffs)
         cffs(1)=ys(j)
         cffs(2)=dys(j)
          do k=1,m
           cffs(k+2)=coeffs(k)
          end do
        xs=x(j+1)
        ys(j+1)=fdf(xx(1),cffs,m2,dys(j+1))
       end do

      j=nmax+1
1     continue
         jlast=j
         id=n-jlast
         m=id-1
           if(m.eq.0) then
            j=j+1
            xs=x(j)
            ys(j)=fdf(xx(2),cffs,m2,dys(j))
            return
           end if
         m2=m+2
         do i=1,id
          xx(i)=x(jlast+i)-xs
          yy(i)=y(jlast+i)-ys(jlast)-dys(jlast)*xx(i)
         end do
         call approx(xx,yy,id,polin2,m,coeffs)
         cffs(1)=ys(jlast)
         cffs(2)=dys(jlast)
          do k=1,m
           cffs(k+2)=coeffs(k)
          end do
        j=j+1
        xs=x(j)
        ys(j)=fdf(xx(1),cffs,m2,dys(j))
      go to 1
      end

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
