! Computation of Lagrangian operators, including:
! time evolution         -> advection()
! Periodic reboxing      -> rebox()
! Linear interpolation   -> interpolation()

! Advection of Lagrangian tracer particles saving in the first field
subroutine advection(xp,dxp,xpold,dxpold,step_dt)
use paran
use commphy
use modlag
implicit none

real(8), intent(inout ), dimension(2,NP) :: xp,xpold
real(8), intent(inout ), dimension(NV,2,NP) :: dxp,dxpold
real(8) :: rebox,step_dt
integer :: ip

! move particle
!$acc kernels present(xp,xpold,step_dt)
xp=xpold+step_dt*xp
!$acc end kernels

! evolution of tangent vectors
!$acc kernels present(dxp,dxpold,step_dt)
dxp=dxpold+step_dt*dxp
!$acc end kernels

! putback particle in the box
!$acc parallel loop collapse(1) present(xp)
do ip=1,NP
 xp(1,ip)=rebox(xp(1,ip),xlx)
 xp(2,ip)=rebox(xp(2,ip),xly)
end do
!$acc end parallel loop

return
end

! Advection of Lagrangian tracer particles saving in second field
subroutine advection2(xp,dxp,xpold,dxpold,step_dt)
use paran
use commphy
use modlag
implicit none

real(8), intent(inout ), dimension(2,NP) :: xp,xpold
real(8), intent(inout ), dimension(NV,2,NP) :: dxp,dxpold
real(8) :: rebox,step_dt
integer :: ip

! move particle
!$acc kernels present(xp,xpold,step_dt)
xpold=xpold+step_dt*xp
!$acc end kernels

! evolution of tangent vectors
!$acc kernels present(dxp,dxpold,step_dt)
dxpold=dxpold+step_dt*dxp
!$acc end kernels

! putback particle in the box
!$acc parallel loop collapse(1) present(xpold)
do ip=1,NP
 xpold(1,ip)=rebox(xpold(1,ip),xlx)
 xpold(2,ip)=rebox(xpold(2,ip),xly)
end do
!$acc end parallel loop

return
end

! --------------------------------------------------------------
! Reboxing function
function rebox(x,xl)
implicit none
real(8) :: x,xl,rebox

rebox=x
if (rebox.lt.0.d0) then
 rebox=rebox+xl
end if

if (rebox.ge.xl) then
 rebox=rebox-xl
end if

return
end 

! --------------------------------------------------------------
! n-th order interpolation (only up to 2 in this version)

subroutine interpolate(u,v,ux,uy,vx,xp,dp)
use paran
use commphy
use modlag
use commsim
implicit none
 
real(8), dimension(NXP2,NY) :: u,v,ux,uy,vx
real(8), intent(inout), dimension(2,NP) :: xp
real(8), intent(inout), dimension(NV,2,NP) :: dp
real(8) :: ci(0:4),cj(0:4)
real(8) :: ul,vl,uxl,uyl,vxl
real(8) :: dx1,dx2
real(8) :: dx,dy,scra
integer i,j,l,ix,iy,ip,iv,il
integer sg,mx,my

!$acc kernels present(u,v,ux,uy,vx,xp,dp,co)
do ip=1,NP
 !temporary variables
 ci=0.d0
 cj=0.d0
 ul=0.d0
 vl=0.d0
 uxl=0.d0
 uyl=0.d0
 vxl=0.d0
 ix=scalax*xp(1,ip) ! Integer of dx in xp
 iy=scalay*xp(2,ip)
 ! Distance for the position and the gridpoint in units of dx
 dx=scalax*xp(1,ip)-real(ix)
 dy=scalay*xp(2,ip)-real(iy) 
 
  ! Interpolation coefficients for the first order (linear)
 ci(0)=(1.0-dx)
 ci(1)=dx
 cj(0)=(1.0-dy)
 cj(1)=dy
 
 !it would be good to save the velocity gradients in the position of the particle
 
! Need to be tested
!  !  Interpolation coefficients for the second order
!  ci(0)=0.5*(1.0-dx)*(2.0-dx)
!  ci(1)=dx*(2.0-dx)
!  ci(2)=-0.5*dx*(1.0-dx)
!  cj(0)=0.5*(1.0-dy)*(2.0-dy)
!  cj(1)=dy*(2.0-dy)
!  cj(2)=-0.5*dy*(1.0-dy)
 
! general interpolation of order oitrp
!  do l=0,oitrp
!  do il=1,oitrp+1
!   ci(l) = ci(l) + co(oitrp,l+1,il) * dx**dble(il-1)
!   cj(l) = cj(l) + co(oitrp,l+1,il) * dy**dble(il-1)
!  enddo
!  enddo
 
 ! Interpolation
 do j=0,oitrp
  do i=0,oitrp
   !find closest square in the grid
   my=mod(iy+j,NY)+1
   mx=mod(ix+i,NX)+1
   scra=ci(i)*cj(j) 
   ! Velocity interpolation
   ul=ul+scra*u(mx,my)
   vl=vl+scra*v(mx,my)
   ! Vel gradient interpolation
   uxl=uxl+scra*ux(mx,my)
   uyl=uyl+scra*uy(mx,my)
   vxl=vxl+scra*vx(mx,my)
   !vyl=-uxl
  enddo
 enddo

! Rewrite velocity
 xp(1,ip)=ul
 xp(2,ip)=vl

 ! Rewrite velocity gradients
 do iv=1,NV
  dx1=dp(iv,1,ip)
  dx2=dp(iv,2,ip)
  dp(iv,1,ip)=uxl*dx1+uyl*dx2
  dp(iv,2,ip)=vxl*dx1-uxl*dx2
 enddo
end do
!$acc end kernels

return
end

