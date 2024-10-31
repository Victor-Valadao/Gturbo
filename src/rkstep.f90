! 2-nd order Runge Kutta with explicit integration of linear terms
subroutine step(z,znl,u,v)
use paran
use commk
use commphy
use commsim
implicit none
complex(8), dimension(NX2P1,NY) :: z,znl,u,v
integer :: i,j

!$acc kernels present(z,znl) ! Copy z to znl
znl=z
!$acc end kernels

! Calculate the nonlinear part
call nlt(znl,u,v) 

! Half step of the time evolution
if(split.eq.1) then
	
	!$acc parallel loop collapse(2) present(z,znl,emkx,emky,dt2)
	do j=1,NY
	do i=1,NX2P1
	 znl(i,j)=emkx(i)*emky(j)*( z(i,j) + dt2*znl(i,j) )
	end do
	end do
	!$acc end parallel 
	
	else
	
	!$acc parallel loop collapse(2) present(z,znl,emk,dt2)
	do j=1,NY
	do i=1,NX2P1
	 znl(i,j)=emk(i,j)*( z(i,j) + dt2*znl(i,j) )
	end do
	end do
	!$acc end parallel 
	
end if

!$acc kernels present(znl)
znl(1,1)=0.d0
!$acc end kernels

! Calculate nonlinear part in half step
call nlt(znl,u,v)

! Last step of the time evolution  
if(split.eq.1) then
	
	!$acc parallel loop collapse(2) present(z,znl,emkx,emky,dt)
	do j=1,NY
	do i=1,NX2P1
	 z(i,j)=emkx(i)*emky(j)*( emkx(i)*emky(j)*z(i,j) + dt*znl(i,j) )
	end do
	end do
	!$acc end parallel 
	
	else
	
	!$acc parallel loop collapse(2) present(z,znl,emk,dt2)
	do j=1,NY
	do i=1,NX2P1
	 z(i,j)=emk(i,j)*( emk(i,j)*z(i,j) + dt*znl(i,j) )
	end do
	end do
	!$acc end parallel 
	
end if

!$acc kernels present(z)
z(1,1)=0.d0
!$acc end kernels

return
end

!======================================================================
! 4-th order Runge Kutta with explicit integration of linear terms


subroutine step4(z,znl,znew,u,v)
use paran
use commk
use commphy
use commsim
implicit none
complex(8), dimension(NX2P1,NY) :: z,znl,znew,u,v

! Step 1
!$acc kernels present(z,znl)
znl=z
!$acc end kernels
call nlt(znl,u,v)
!$acc kernels present(z,znl,znew,emk,dt6)
znew=emk*emk*(z+dt6*znl)
!$acc end kernels

! Step 2
!$acc kernels present(z,znl,emk,dt2)
znl=emk*(z+dt2*znl)
!$acc end kernels
call nlt(znl,u,v)
!$acc kernels present(znl,znew,emk,dt3)
znew=znew+dt3*emk*znl
!$acc end kernels

! Step 3
!$acc kernels present(z,znl,emk,dt2)
znl=emk*z+dt2*znl
!$acc end kernels
call nlt(znl,u,v)
!$acc kernels present(znl,znew,emk,dt3)
znew=znew+dt3*emk*znl
!$acc end kernels

! Step 4
!$acc kernels present(z,znl,emk,dt)
znl=emk*(emk*z+dt*znl)
!$acc end kernels
call nlt(znl,u,v)
!$acc kernels present(z,znl,znew,dt6)
z=znew+dt6*znl
!$acc end kernels

!$acc kernels present(z)
z(1,1)=0.d0
!$acc end kernels

return
end


subroutine step1(z,znl)
use paran
use commk
use commphy
use commsim
implicit none

complex(8), intent(in    ), dimension(NX2P1,NY) :: z
complex(8), intent(inout ), dimension(NX2P1,NY) :: znl
integer i,j


!$acc parallel loop collapse(2) present(z,znl,emk,dt2)
do j=1,NY
do i=1,NX2P1
 znl(i,j)=emk(i,j)*(z(i,j)+dt2*znl(i,j))
end do
end do
!$acc end parallel 

return
end

subroutine step2(z,znl)
use paran
use commk
use commphy
use commsim
implicit none

complex(8), intent(inout ), dimension(NX2P1,NY) :: z
complex(8), intent(in    ), dimension(NX2P1,NY) :: znl
integer :: i,j

!$acc parallel loop collapse(2) present(z,znl,emk,dt)
do j=1,NY
do i=1,NX2P1
 z(i,j)=emk(i,j)*(emk(i,j)*z(i,j)+dt*znl(i,j))
end do
end do
!$acc end parallel

return
end
