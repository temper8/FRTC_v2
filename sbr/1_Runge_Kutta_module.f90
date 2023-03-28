module Runge_Kutta_module
    implicit none
    
contains
    
!----------------------------------------------------------------

    subroutine rkqc(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      implicit real*8 (a-h,o-z)
      integer, intent(in) :: n
      integer nmax, i
      parameter (nmax=10,fcor=.0666666667d0, &
         one=1.d0,safety=0.9d0,errcon=6.d-4)
      external derivs
      dimension y(n),dydx(n),yscal(n),ytemp(nmax),ysav(nmax),dysav(nmax)
      pgrow=-0.20d0
      pshrnk=-0.25d0
      xsav=x
      do 11 i=1,n
        ysav(i)=y(i)
        dysav(i)=dydx(i)
11    continue
      h=htry
1     hh=0.5d0*h
      call rk4(ysav,dysav,n,xsav,hh,ytemp,derivs)
      x=xsav+hh
      call derivs(x,ytemp,dydx)
      call rk4(ytemp,dydx,n,x,hh,y,derivs)
      x=xsav+h
      if(x.eq.xsav) then
       write(*,*)' stepsize not significant in rkqc'
       write(*,*)'xsav=',xsav,' h=',h,' htry=',htry
       write(*,88) y,dydx
       pause
      end if
88    format(1x,10(e14.7,1x))

      call rk4(ysav,dysav,n,xsav,h,ytemp,derivs)
      errmax=0.d0
      do 12 i=1,n
        v=ytemp(i)
        ytemp(i)=y(i)-ytemp(i)
        errmax=dmax1(errmax,dabs(ytemp(i)/yscal(i)))
12    continue
      errmax=errmax/eps
      if(errmax.gt.one) then
        h=safety*h*(errmax**pshrnk)
        goto 1
      else
        hdid=h
        if(errmax.gt.errcon)then
          hnext=safety*h*(errmax**pgrow)
        else
          hnext=4.d0*h
        endif
      endif
      do 13 i=1,n
        y(i)=y(i)+ytemp(i)*fcor
13    continue
      return
     end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rk4(y,dydx,n,x,h,yout,derivs)
      implicit real*8 (a-h,o-z)
      integer, intent(in) :: n
      integer nmax, i
      parameter (nmax=10)
      dimension y(n),dydx(n),yout(n),yt(nmax),dyt(nmax),dym(nmax)
      hh=h*0.5d0
      h6=h/6.d0
      xh=x+hh
      do 11 i=1,n
        yt(i)=y(i)+hh*dydx(i)
11    continue
      call derivs(xh,yt,dyt)
        dv1=dyt(3)
      do 12 i=1,n
        yt(i)=y(i)+hh*dyt(i)
12    continue
      call derivs(xh,yt,dym)
      do 13 i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
13    continue

      call derivs(x+h,yt,dyt)
      do 14 i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.d0*dym(i))
14    continue
      return
    end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!sav2008: below this line there are new subroutins and functions
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
        integer n,nmax
        double precision eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
        external derivs
        parameter (nmax=50)
  !cu    uses derivs,rkck
        integer i
        double precision errmax,h,htemp,xnew,yerr(nmax),ytemp(nmax),safety,pgrow,pshrnk,errcon
        parameter (safety=0.9d0,pgrow=-.2d0,pshrnk=-.25d0,errcon=1.89d-4)
        h=htry
  1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
        errmax=0.d0
        do i=1,n
          errmax=max(errmax,abs(yerr(i)/yscal(i)))
        enddo
        errmax=errmax/eps
        if(errmax.gt.1.d0)then
          htemp=safety*h*(errmax**pshrnk)
          h=sign(max(abs(htemp),0.1d0*abs(h)),h)
          xnew=x+h
          if(xnew.eq.x)pause 'stepsize underflow in rkqs'
          goto 1
        else
          if(errmax.gt.errcon)then
            hnext=safety*h*(errmax**pgrow)
          else
            hnext=5.d0*h
          endif
          hdid=h
          x=x+h
          do i=1,n
            y(i)=ytemp(i)
          enddo
          return
        endif
    end
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rkck(y,dydx,n,x,h,yout,yerr,derivs)
        integer n,nmax
        double precision h,x,dydx(n),y(n),yerr(n),yout(n)
        external derivs
        parameter (nmax=50)
  !cu    uses derivs
        integer i
        double precision ak2(nmax),ak3(nmax),ak4(nmax),ak5(nmax),ak6(nmax) &
            ,ytemp(nmax),a2,a3,a4,a5,a6,b21,b31,b32,b41,b42,b43,b51,b52,b53, &
            b54,b61,b62,b63,b64,b65,c1,c3,c4,c6,dc1,dc3,dc4,dc5,dc6
        parameter (a2=.2d0,a3=.3d0,a4=.6d0, a5=1.d0, a6=.875d0, b21=.2d0, b31=3.d0/40.d0, &
            b32=9.d0/40.d0,b41=.3d0,b42=-.9d0,b43=1.2d0,b51=-11.d0/54.d0, b52=2.5d0, &
            b53=-70.d0/27.d0,b54=35.d0/27.d0,b61=1631.d0/55296.d0, b62=175.d0/512.d0, &
            b63=575.d0/13824.d0,b64=44275.d0/110592.d0,b65=253.d0/4096.d0, c1=37.d0/378.d0, &
            c3=250.d0/621.d0,c4=125.d0/594.d0,c6=512.d0/1771.d0,dc1=c1-2825.d0/27648.d0, &
            dc3=c3-18575.d0/48384.d0,dc4=c4-13525.d0/55296.d0,dc5=-277.d0/14336.d0, dc6=c6-.25d0)
        do i=1,n
          ytemp(i)=y(i)+b21*h*dydx(i)
        enddo
        call derivs(x+a2*h,ytemp,ak2)
        do i=1,n
          ytemp(i)=y(i)+h*(b31*dydx(i)+b32*ak2(i))
        enddo
        call derivs(x+a3*h,ytemp,ak3)
        do i=1,n
          ytemp(i)=y(i)+h*(b41*dydx(i)+b42*ak2(i)+b43*ak3(i))
        enddo
        call derivs(x+a4*h,ytemp,ak4)
        do i=1,n
          ytemp(i)=y(i)+h*(b51*dydx(i)+b52*ak2(i)+b53*ak3(i)+b54*ak4(i))
        enddo
        call derivs(x+a5*h,ytemp,ak5)
        do  i=1,n
          ytemp(i)=y(i)+h*(b61*dydx(i)+b62*ak2(i)+b63*ak3(i)+b64*ak4(i)+b65*ak5(i))
        enddo
        call derivs(x+a6*h,ytemp,ak6)
        do  i=1,n
          yout(i)=y(i)+h*(c1*dydx(i)+c3*ak3(i)+c4*ak4(i)+c6*ak6(i))
        enddo
        do  i=1,n
          yerr(i)=h*(dc1*dydx(i)+dc3*ak3(i)+dc4*ak4(i)+dc5*ak5(i)+dc6*ak6(i))
        enddo
        return
        end
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
end module Runge_Kutta_module