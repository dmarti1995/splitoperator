PROGRAM EULERS
integer,parameter:: N=1200,N2=700
real(8):: dt,time,tf,hbar,m,ome0,pi,t0,k0
real(8):: xr(N),xi(N),dx,x(N),V(N),L,alpha
real(8):: sigma,x0,t1,fr(N),fi(N),kf,k00,hk
real(8):: frw(N2+1),fiw(N2+1),k,norma,normak,ks(N2+1)
integer:: ii,nts,jj,iter,kkk,jjj
character(len=13):: filename
hbar=1.0
x0=0.0
pi=acos(-1.0)
m=1.0
V=0.0
L=20.0
dx=L/(N-1)
print *, dx
open(11,file='potential.dat')
do ii=1,N
   x(ii)=-10+(ii-1)*dx
!   if (x(ii)<4.0) then
      V(ii)=0.0
     ! V(ii)=-48+3*x(ii)**2
  ! else
   !   V(ii)=150.0
   !endif
   write(11,*) x(ii),V(ii)
enddo
close(11)
dt=1D-5
tf=2.0
nts=int(tf/dt)
time=0.0
alpha=1.0
sigma=1.0/alpha
t0=0.0
t1=-dt/2.0
k0=10.0
ome0=k0**2/(2*m)
call initial(N,sigma,pi,x0,xi,xr,x,k0,t0,t1,ome0)
ind=int(N/2)
iter=1
open(7,file='fourier.dat')
kf=k0+10
k00=k0-30
hk=(kf-k00)/N2
do kkk=1,N2+1
   ks(kkk)=k00+(kkk-1)*hk
enddo   
open(8,file='eulerint.dat')
do ii=1,nts
   if (mod(ii,100)==0) then
       call writing(N,N2,x,ks,xr,xi,fiw,frw)
       iter=iter+1
       print *, iter
   endif
   !normak=0.0
   !do kkk=1,1001
   !   k=k00+(kkk-1)*hk
   !   do jjj=1,N
!         fr(jjj)=xr(jjj)*cos(k*x(jjj))+xi(jjj)*sin(k*x(jjj))
!         fi(jjj)=xi(jjj)*cos(k*x(jjj))-xr(jjj)*sin(k*x(jjj))
!      enddo
!      call simpson(fr,dx,N,frw(kkk))
!      call simpson(fi,dx,N,fiw(kkk))
!      fiw(kkk)=fiw(kkk)/sqrt(2*pi)
!      frw(kkk)=frw(kkk)/sqrt(2*pi)
!      write(7,*) k, frw(kkk), fiw(kkk), frw(kkk)**2+fiw(kkk)**2
!      normak=normak+(frw(kkk)**2+fiw(kkk)**2)*hk
!    enddo        
!      write(7,*)
!      write(7,*)
!      print *, 'k',normak
!   norma=0.0
!   do jj=1,N
!      norma=norma+dx*(xi(jj)**2+xr(jj)**2)
!      write(8,*) x(jj),xr(jj),xi(jj),xi(jj)**2+xr(jj)**2
!   enddo
!   write(8,*)
!   write(8,*)
!   print *, 'x',norma
!   endif
!   CALL rk4(time,xr,xi,N,dt,dx,V,hbar,m)
   CALL splitoperator(N,N2,xr,xi,fiw,frw,dx,time,ks,hk,x,V,m)
   time=time+dt
enddo
   do kkk=1,N2+1
      k=k00+(kkk-1)*hk
      do jjj=1,N
         fr(jjj)=xr(jjj)*cos(k*x(jjj))+xi(jjj)*sin(k*x(jjj))
         fi(jjj)=xi(jjj)*cos(k*x(jjj))-xr(jjj)*sin(k*x(jjj))
      enddo
      call simpson(fr,dx,N,frw(kkk))
      call simpson(fi,dx,N,fiw(kkk))
      fiw(kkk)=fiw(kkk)/sqrt(2*pi)
      frw(kkk)=frw(kkk)/sqrt(2*pi)
      write(7,*) k, frw(kkk), fiw(kkk), frw(kkk)**2+fiw(kkk)**2
    enddo        
do jj=1,N
   write(8,*) x(jj),xr(jj),xi(jj),xr(jj)**2+xi(jj)**2
enddo
close(8)
close(7)

contains

subroutine initial(N,sigma,pi,x0,xi,xr,x,k0,t0,t1,ome0)
integer:: N,ii
real(8):: sigma,pi,x0,k0,t0,ome0,t1
real(8):: xi(N),xr(N),x(N)
do ii=1,N
   a=(1.0/sqrt(sigma*sqrt(pi)))*exp((-(x(ii)-x0)**2)/(2.0*sigma**2))
   xr(ii)=a*cos(k0*x(ii)-ome0*t0)
   xi(ii)=a*sin(k0*x(ii)-ome0*t1)
enddo
xr(1)=0.0
xi(1)=0.0
xr(N)=0.0
xi(N)=0.0
endsubroutine

subroutine FTCS(ti,xr,xi,N,dt,dx,V,hbar,m)
integer:: N
real(8):: xr(N),ti,dt,xr2(N),xi(N),xi2(N),V(N)
real(8):: hbar,dx,m
!xi2(1)=xi(1)+(dt*hbar/(dx**2.0*2*m))*(xr(2)-xr(1))-(dt/hbar)*V(1)*xr(1)
!xi2(N)=xi(N)+(dt*hbar/(dx**2.0*2*m))*(xr(N-1)-xr(N))-(dt/hbar)*V(N)*xr(N)
do jj=2,N-1
   xi2(jj)=xi(jj)-dt*(xr(jj)*V(jj)/hbar-hbar*(xr(jj+1)+xr(jj-1)-2*xr(jj))/(2*m*dx**2))
enddo

!xr2(1)=xr(1)-(dt*hbar/(dx**2.0*2*m))*(xi(2)-xi(1))+(dt/hbar)*V(1)*xi(1)
!:qxr2(N)=xr(N)-(dt*hbar/(dx**2.0*2*m))*(xi(N-1)-xi(N))+(dt/hbar)*V(N)*xi(N)
do jj=2,N-1
   xr2(jj)=xr(jj)+dt*(xi(jj)*V(jj)/hbar-hbar*(xi(jj+1)+xi(jj-1)-2*xi(jj))/(2*m*dx**2))
enddo
xi=xi2
xr=xr2
end subroutine

subroutine rk4(ti,xr,xi,N,dt,dx,V,hbar,m)
integer:: N
real(8):: xr(N),ti,dt,xr2(N),xi(N),xi2(N),V(N)
real(8):: hbar,dx,m,xr3(N),xi3(N),xr4(N),xi4(N)
real(8):: k1i(N),k1r(N),k2i(N),k2r(N),k3i(N),k3r(N)
real(8):: k4r(N),k4i(N)
do jj=2,N-1
   k1i(jj)=-(xr(jj)*V(jj)/hbar-hbar*(xr(jj+1)+xr(jj-1)-2*xr(jj))/(2*m*dx**2))
   k1r(jj)=(xi(jj)*V(jj)/hbar-hbar*(xi(jj+1)+xi(jj-1)-2*xi(jj))/(2*m*dx**2))
   xi2(jj)=xi(jj)+(dt*k1i(jj)/2.0)
   xr2(jj)=xr(jj)+(dt*k1r(jj)/2.0)
enddo

do jj=2,N-1
   k2i(jj)=-(xr2(jj)*V(jj)/hbar-hbar*(xr2(jj+1)+xr2(jj-1)-2*xr2(jj))/(2*m*dx**2))
   k2r(jj)=(xi2(jj)*V(jj)/hbar-hbar*(xi2(jj+1)+xi2(jj-1)-2*xi2(jj))/(2*m*dx**2))
   xi3(jj)=xi2(jj)+(dt*k2i(jj)/2.0)
   xr3(jj)=xr2(jj)+(dt*k2r(jj)/2.0)
enddo

do jj=2,N-1
   k3i(jj)=-(xr3(jj)*V(jj)/hbar-hbar*(xr3(jj+1)+xr3(jj-1)-2*xr3(jj))/(2*m*dx**2))
   k3r(jj)=(xi3(jj)*V(jj)/hbar-hbar*(xi3(jj+1)+xi3(jj-1)-2*xi3(jj))/(2*m*dx**2))
   xi4(jj)=xi3(jj)+(dt*k3i(jj))
   xr4(jj)=xr3(jj)+(dt*k3r(jj))
enddo

do jj=2,N-1
   k4i(jj)=-(xr4(jj)*V(jj)/hbar-hbar*(xr4(jj+1)+xr4(jj-1)-2*xr4(jj))/(2*m*dx**2))
   k4r(jj)=(xi4(jj)*V(jj)/hbar-hbar*(xi4(jj+1)+xi4(jj-1)-2*xi4(jj))/(2*m*dx**2))
   xi(jj)=xi(jj)+(dt/6.0)*(k1i(jj)+2.0*k2i(jj)+2.0*k3i(jj)+k4i(jj))
   xr(jj)=xr(jj)+(dt/6.0)*(k1r(jj)+2.0*k2r(jj)+2.0*k3r(jj)+k4r(jj))
enddo
end subroutine

SUBROUTINE simpson(valorsfunc,dist,N,resultat)
integer:: N
real(8),dimension(N):: valorsfunc
real(8):: dist,resultat
resultat=(dist/3)*(valorsfunc(1)+valorsfunc(N))
do i=2,(N-1)
if (mod(i,2)==0) then
resultat=resultat+(4.0*dist/3.0)*valorsfunc(i)
else
resultat=resultat+(2.0*dist/3.0)*valorsfunc(i)
end if
end do
return
END SUBROUTINE

SUBROUTINE splitoperator(N,N2,xr,xi,fiw,frw,dx,dt,ks,hk,x,V,m)
integer:: N,N2
real(8):: xr(N),xi(N),fiw(N2+1),frw(N2+1),dx,dt,k00,hk,x(N),V(N),m
real(8):: fr(N),fi(N),k,fr2(N2+1),fi2(N2+1),ks(N2+1),fiw2(N),frw2(N)
integer::kkk,jjj
do kkk=1,N2+1
   k=ks(kkk)
   do jjj=1,N
      fr(jjj)=cos(k*x(jjj))*(xr(jjj)*cos(V(jjj)*dt/2.0)+xi(jjj)*sin(V(jjj)*dt/2.0))+&
              sin(k*x(jjj))*(-xr(jjj)*sin(V(jjj)*dt/2.0)+xi(jjj)*cos(V(jjj)*dt/2.0))
      fi(jjj)=cos(k*x(jjj))*(-xr(jjj)*sin(V(jjj)*dt/2.0)+xi(jjj)*cos(V(jjj)*dt/2.0))+&
              sin(k*x(jjj))*(-xr(jjj)*cos(V(jjj)*dt/2.0)-xi(jjj)*sin(V(jjj)*dt/2.0))
   enddo
   call simpson(fr,dx,N,frw(kkk))
   call simpson(fi,dx,N,fiw(kkk))
   fiw(kkk)=fiw(kkk)/sqrt(2*pi)
   frw(kkk)=frw(kkk)/sqrt(2*pi)
enddo
do jjj=1,N
   do kkk=1,N2+1
      k=ks(kkk)
      fr2(kkk)=cos(k*x(jjj))*(xr(jjj)*cos(k**2.0*dt/(2.0*m))+xi(jjj)*sin(k**2.0*dt/(2.0*m)))-&
              sin(k*x(jjj))*(-xr(jjj)*sin(k**2.0*dt/(2.0*m))+xi(jjj)*cos(k**2.0*dt/(2.0*m)))
      fi2(kkk)=cos(k*x(jjj))*(-xr(jjj)*sin(k**2.0*dt/(2.0*m))+xi(jjj)*cos(k**2.0*dt/(2.0*m)))-&
              sin(k*x(jjj))*(-xr(jjj)*cos(k**2.0*dt/(2.0*m))-xi(jjj)*sin(k**2.0*dt/(2.0*m)))
   enddo
   call simpson(fr2,hk,N2+1,frw(jjj))
   call simpson(fi2,hk,N2+1,fiw(jjj))
   fiw2(jjj)=fiw(jjj)/sqrt(2*pi)
   frw2(jjj)=frw(jjj)/sqrt(2*pi)
   xr(jjj)=cos(V(jjj)*dt/2.0)*frw2(jjj)+sin(V(jjj)*dt/2.0)*fiw2(jjj)
   xi(jjj)=-sin(V(jjj)*dt/2.0)*frw2(jjj)+cos(V(jjj)*dt/2.0)*fiw2(jjj)
enddo         
END SUBROUTINE

SUBROUTINE writing(N,N2,x,ks,xr,xi,fiw,frw)
integer:: N,N2
real(8):: x(N),xr(N),xi(N),ks(N2+1),fiw(N2+1),frw(N2+1)
do jj=1,N
   write(8,*) x(jj),xr(jj),xi(jj),xi(jj)**2+xr(jj)**2
enddo
write(8,*)
write(8,*)
do kkk=1,N2+1
   write(7,*) ks(kkk), frw(kkk), fiw(kkk), frw(kkk)**2+fiw(kkk)**2
enddo        
write(7,*)
write(7,*)
END SUBROUTINE

END PROGRAM
