! Some lyap functions, including:
! rescale()   -> rescaling of the pert field
! add_noise() -> initial perturbation

subroutine rescale(z,w,eed)

use paran
use commk
use lyap

complex(8), dimension(NX2P1,NY) :: z,w
real*8 eed,k2,scra,fac
integer i,j,k,ierr


!$acc kernels present(eed) !
eed=0.d0
!$acc end kernels

! dw(t)=(w(t)-z(t))/sqrt(2)
!$acc parallel loop collapse(2) present(z,w,eed)
do j=1,NY
do i=1,NX2P1
 fac=1.0
 if ((i.eq.1).or.(i.eq.NX2P1)) fac=0.5d0
 scra=0.5d0*cdabs(z(i,j)-w(i,j))**2 
 eed=eed+fac*scra      ! Z[dw](t)
end do
end do
!$acc end parallel loop

!$acc kernels present(eed) !
eed=delta0/sqrt(2.d0*eed)  ! factor=delta0/sqrt(2Z[dw](t))
!$acc end kernels

! w'(0) = z(t)+factor*(w(t)-z(t))
! dw'(0)= (w'(0)-z(t))/sqrt(2)
! Z[dw'](0) = factor**2 Z[dw](t) = delta0**2/2
! Z[dw'](0) = Z[dw](0)
!$acc parallel loop collapse(2) present(z,w,eed)
do j=1,NY
do i=1,NX2P1
 w(i,j)=z(i,j)+eed*(w(i,j)-z(i,j)) 
end do
end do
!$acc end parallel loop

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
integer :: i,j,seed,i0,j0
real :: rann,a2

! a = 2 sqrt(2) sqrt(3): corr, def dw, trunc
! white space-time noise

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

end subroutine
