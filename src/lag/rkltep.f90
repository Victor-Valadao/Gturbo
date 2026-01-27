! Steps on time evolution, including:
! step()   <- Rk2
! step4()  <- Rk4
! Split of Emk matrix deprecated in this version
! USE OF ADVECTION OR ADVECTION2 DEPENDS ON WHERE YOU SAVE THE RESULT

!=================================================================
! 2-nd order Runge Kutta with explicit integration of linear terms
subroutine step(z,znl,u,v,dux,duy,dvx,xp,dxp,xpnl,dxpnl)
use paran
use commk
use commphy
use commsim
use modlag
implicit none
complex(8), dimension(NX2P1,NY) :: z,znl,u,v,dux,duy,dvx
real(8), dimension(NV,2,NP) :: dxp,dxpnl
real(8), dimension(2,NP) :: xp,xpnl

!$acc kernels present(z,znl) ! Copy z to znl
znl=z
!$acc end kernels
!$acc kernels present(xp,xpnl)
xpnl=xp
!$acc end kernels
!$acc kernels present(dxp,dxpnl)
dxpnl=dxp
!$acc end kernels

! Calculate the nonlinear part 
call nlt(znl,u,v,dux,duy,dvx,xpnl,dxpnl)
call advection(xpnl,dxpnl,xp,dxp,dt2)

!$acc kernels present(z,znl,emk,dt2)
znl=emk*(z + dt2*znl)
!$acc end kernels

! Calculate nonlinear part in half step
call nlt(znl,u,v,dux,duy,dvx,xpnl,dxpnl)
call advection(xpnl,dxpnl,xp,dxp,dt)

!$acc kernels present(xp,xpnl)
xp=xpnl
!$acc end kernels

!$acc kernels present(dxp,dxpnl)
dxp=dxpnl
!$acc end kernels

!$acc kernels present(z,znl,emk,dt)
z=emk*(emk*z + dt*znl)
!$acc end kernels

!$acc kernels present(z)
z(1,1)=0.d0
!$acc end kernels

return
end

!======================================================================
! 4-th order Runge Kutta with explicit integration of linear terms

subroutine step4(z,znl,znew,u,v,dux,duy,dvx,xp,dxp,xpnl,dxpnl,xpnew,dxpnew)
use paran
use commk
use commphy
use commsim
use modlag
implicit none
complex(8), dimension(NX2P1,NY) :: z,znl,znew,u,v,dux,duy,dvx
real(8), dimension(NV,2,NP) :: dxp,dxpnl,dxpnew
real(8), dimension(2,NP) :: xp,xpnl,xpnew

! Step 1
!$acc kernels present(z,znl) ! Copy z to znl
znl=z
!$acc end kernels
!$acc kernels present(xp,xpnl,xpnew)
xpnl=xp
xpnew=xp
!$acc end kernels
!$acc kernels present(dxp,dxpnl,dxpnew)
dxpnl=dxp
dxpnew=dxp
!$acc end kernels
call nlt(znl,u,v,dux,duy,dvx,xpnl,dxpnl)
call advection2(xpnl,dxpnl,xpnew,dxpnew,dt6)
!$acc kernels present(z,znl,znew,emk,dt6)
znew=emk*emk*(z+dt6*znl)
!$acc end kernels

! Step 2
!$acc kernels present(z,znl,emk,dt2)
znl=emk*(z+dt2*znl)
!$acc end kernels
!$acc kernels present(xp,xpnl,dt2)
xpnl=xp+dt2*xpnl
!$acc end kernels
!$acc kernels present(dxp,dxpnl,dt2)
dxpnl=dxp+dt2*dxpnl
!$acc end kernels
call nlt(znl,u,v,dux,duy,dvx,xpnl,dxpnl)
call advection2(xpnl,dxpnl,xpnew,dxpnew,dt3)
!$acc kernels present(znl,znew,emk,dt3)
znew=znew+dt3*emk*znl
!$acc end kernels

! Step 3
!$acc kernels present(z,znl,emk,dt2)
znl=emk*z+dt2*znl
!$acc end kernels
!$acc kernels present(xp,xpnl,dt2)
xpnl=xp+dt2*xpnl
!$acc end kernels
!$acc kernels present(dxp,dxpnl,dt2)
dxpnl=dxp+dt2*dxpnl
!$acc end kernels
call nlt(znl,u,v,dux,duy,dvx,xpnl,dxpnl)
call advection2(xpnl,dxpnl,xpnew,dxpnew,dt3)
!$acc kernels present(znl,znew,emk,dt3)
znew=znew+dt3*emk*znl
!$acc end kernels

! Step 4
!$acc kernels present(z,znl,emk,dt)
znl=emk*(emk*z+dt*znl)
!$acc end kernels
!$acc kernels present(xp,xpnl,dt)
xpnl=xp+dt*xpnl
!$acc end kernels
!$acc kernels present(dxp,dxpnl,dt)
dxpnl=dxp+dt*dxpnl
!$acc end kernels
call nlt(znl,u,v,dux,duy,dvx,xpnl,dxpnl)
call advection2(xpnl,dxpnl,xpnew,dxpnew,dt6)
!$acc kernels present(z,znl,znew,dt6)
z=znew+dt6*znl
!$acc end kernels

!$acc kernels present(xp,xpnew)
xp=xpnew
!$acc end kernels

!$acc kernels present(dxp,dxpnew)
dxp=dxpnew
!$acc end kernels

!$acc kernels present(z)
z(1,1)=0.d0
!$acc end kernels

return
end