! Steps on time evolution, including:
! step()   <- Rk2
! step4()  <- Rk4
! Split of Epk matrix deprecated in this version

!version for the passive scalar with diff params

!======================================================================
! 2-nd order Runge Kutta with explicit integration of linear terms
subroutine step(z,znl,w,wnl,u,v)
use paran
use commk
use commphy
use commsim
implicit none
complex(8), dimension(NX2P1,NY) :: z,znl
complex(8), dimension(NX2P1,NY) :: w,wnl
complex(8), dimension(NX2P1,NY) :: u,v
integer :: i,j

!$acc kernels present(z,w,znl,wnl) ! Copy z to znl
znl=z
wnl=w
!$acc end kernels

! Calculate the nonlinear part
call nlt(znl,wnl,u,v)

! Half step of the time evolution
!$acc parallel loop collapse(2) present(z,znl,emk,w,wnl,epk,dt2)
do j=1,NY
do i=1,NX2P1
 znl(i,j)=emk(i,j)*( z(i,j) + dt2*znl(i,j) )
 wnl(i,j)=epk(i,j)*( w(i,j) + dt2*wnl(i,j) )
end do
end do
!$acc end parallel 

!$acc kernels present(znl,wnl)
znl(1,1)=0.d0
wnl(1,1)=0.d0
!$acc end kernels

! Calculate nonlinear part in half step
call nlt(znl,wnl,u,v)
	
!$acc parallel loop collapse(2) present(z,znl,emk,w,wnl,epk,dt)
do j=1,NY
do i=1,NX2P1
 z(i,j)=emk(i,j)*( emk(i,j)*z(i,j) + dt*znl(i,j) )
 w(i,j)=epk(i,j)*( epk(i,j)*w(i,j) + dt*wnl(i,j) )
end do
end do
!$acc end parallel 
	
!$acc kernels present(z,w)
z(1,1)=0.d0
w(1,1)=0.d0
!$acc end kernels

return
end

!======================================================================
! 4-th order Runge Kutta with explicit integration of linear terms

subroutine step4(z,znl,znew,w,wnl,wnew,u,v)
use paran
use commk
use commphy
use commsim
implicit none
complex(8), dimension(NX2P1,NY) :: z,znl,znew
complex(8), dimension(NX2P1,NY) :: w,wnl,wnew
complex(8), dimension(NX2P1,NY) :: u,v

! Step 1
!$acc kernels present(z,znl,w,wnl)
znl=z
wnl=w
!$acc end kernels
call nlt(znl,wnl,u,v)
!$acc kernels present(z,znl,znew,emk,w,wnl,wnew,epk,dt6)
znew=emk*emk*( z + dt6*znl)
wnew=epk*epk*( w + dt6*wnl)
!$acc end kernels

! Step 2
!$acc kernels present(z,znl,emk,w,wnl,epk,dt2)
znl=emk*( z + dt2*znl)
wnl=epk*( w + dt2*wnl)
!$acc end kernels
call nlt(znl,wnl,u,v)
!$acc kernels present(znl,znew,emk,wnl,wnew,epk,dt3)
znew=znew + dt3*emk*znl
wnew=wnew + dt3*epk*wnl
!$acc end kernels

! Step 3
!$acc kernels present(z,znl,emk,w,wnl,epk,dt2)
znl=emk*z + dt2*znl
wnl=epk*w + dt2*wnl
!$acc end kernels
call nlt(znl,wnl,u,v)
!$acc kernels present(znl,znew,emk,wnl,wnew,epk,dt3)
znew=znew + dt3*emk*znl
wnew=wnew + dt3*epk*wnl
!$acc end kernels

! Step 4
!$acc kernels present(z,znl,emk,w,wnl,epk,dt)
znl=emk*(emk*z + dt*znl)
wnl=epk*(epk*w + dt*wnl)
!$acc end kernels
call nlt(znl,wnl,u,v)
!$acc kernels present(z,znl,znew,w,wnl,wnew,dt6)
z=znew + dt6*znl
w=wnew + dt6*wnl
!$acc end kernels

!$acc kernels present(z,w)
z(1,1)=0.d0
w(1,1)=0.d0
!$acc end kernels

return
end
