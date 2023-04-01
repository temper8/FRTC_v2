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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
