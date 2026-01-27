! Some lyap functions, including:
! rescale()   -> rescaling of the perturbed field
! add_noise() -> eulerian perturbation (why?)


! --------------------------------------------------------------
! Gram-Schmidt orthonormalization
subroutine rescale(dxp)
use modlag
use commsim
use commphy
implicit none

real(8) :: norm,prod,scra
real(8), intent(inout ), dimension(NV,2,NP) :: dxp
integer :: ip,i,j,ipas

!$acc kernels present(lyap,ftle)
lyap=0.d0
ftle=0.d0
!$acc end kernels

norm=1.d0/dble(NP)
!$acc kernels present(dxp,lyap,ftle)
do ip=1,NP
  ! max lyap
  scra = dxp(1,1,ip)*dxp(1,1,ip) + dxp(1,2,ip)*dxp(1,2,ip)
  
  ftle(ip,1)=ftle(ip,1)+dlog(scra)/2.d0    ! each particle
  lyap(1)=lyap(1)+norm*dlog(scra)/2.d0     ! mean over the particles
  
  ! min lyap
  prod = dxp(1,1,ip)*dxp(2,1,ip) + dxp(1,2,ip)*dxp(2,2,ip)
  
  ! orthogonal: scra = |u|^2; prod = u.v
  dxp(2,1,ip) = dxp(2,1,ip) - prod*dxp(1,1,ip)/scra
  dxp(2,2,ip) = dxp(2,2,ip) - prod*dxp(1,2,ip)/scra
  
  ! prod -> |v'|^2
  prod = dxp(2,1,ip)*dxp(2,1,ip) + dxp(2,2,ip)*dxp(2,2,ip)
  
  ftle(ip,2)=ftle(ip,2)+dlog(prod)/2.d0    ! each particle
  lyap(2)=lyap(2)+norm*dlog(prod)/2.d0     ! mean over the particles
	
  ! normalization to norm L2 = 1 
  dxp(1,1,ip) = dxp(1,1,ip)/dsqrt(scra)
  dxp(1,2,ip) = dxp(1,2,ip)/dsqrt(scra)
  dxp(2,1,ip) = dxp(2,1,ip)/dsqrt(prod)
  dxp(2,2,ip) = dxp(2,2,ip)/dsqrt(prod)

end do
!$acc end kernels

return
end

! --------------------------------------------------------------
! Gram-Schmidt orthonormalization  !! OLD OUTSIDE GPU
subroutine rescale2(dxp)
use modlag
use commsim
use commphy
implicit none

real(8) :: norm,prod,scra
real(8), intent(inout ), dimension(NV,2,NP) :: dxp
integer :: ip,i,j,ipas

lyap=0.d0
ftle=0.d0

norm=1.d0/dble(NP)
! $acc parallel loop collapse(2) present(dxp,lyap) reduction(+:lyap) reduction(+:ftle)
do ip=1,NP
do i=1,nv

 do j=1,i-1
  !G.S.
  prod = dxp(i,1,ip)*dxp(j,1,ip) + dxp(i,2,ip)*dxp(j,2,ip)
  dxp(i,1,ip) = dxp(i,1,ip) - prod*dxp(j,1,ip)
  dxp(i,2,ip) = dxp(i,2,ip) - prod*dxp(j,2,ip)
 end do
 
 ! |Dx|
 prod=dsqrt( dxp(i,1,ip)*dxp(i,1,ip) + dxp(i,2,ip)*dxp(i,2,ip) )
 
 ftle(ip,i)=ftle(ip,i)+dlog(prod)
 lyap(i)=lyap(i)+norm*dlog(prod)
 
 scra=1.d0/prod
 dxp(i,1,ip)=dxp(i,1,ip)*scra
 dxp(i,2,ip)=dxp(i,2,ip)*scra

end do
end do
! $acc end parallel loop

return
end

!=================================================================

subroutine add_noise(w,seed)
use paran
use lyap
use commphy
implicit none
real(8), intent(inout), dimension(NXP2,NY) :: w
real(8), parameter :: a = sqrt(24.0d0)
integer :: i,j,seed
real :: rann

! a = 2 sqrt(2) sqrt(3): corr, def dw, trunc

do j = 1, NY
do i = 1, NX
w(i,j) = delta0*a*(rann(seed)-0.5d0)
end do
end do

do j=1,NY
do i=NX+1,NXP2
 w(i,j)=0.d0
end do
end do

! call writep(w,1)

end subroutine



