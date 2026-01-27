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

a2=2!*(6.2831853071795864769d0/NX)
! i0=int(rann(seed)*NX/2)
! j0=int(rann(seed)*NY/2)

i0=NX/4
j0=NY/4

do j = 1, NY
do i = 1, NX
w(i,j) = delta0*( dexp(-0.5d0*( (i-i0)**2 + (j-j0)**2 )/(a2*a2)) )
end do
end do

write(6,*)i0,j0
call flush(6)

! do j = 1, NY
! do i = 1, NX
! w(i,j) = delta0*a*(rann(seed)-0.5d0)
! end do
! end do

do j=1,NY
do i=NX+1,NXP2
 w(i,j)=0.d0
end do
end do

! call writep(w,1)

end subroutine


subroutine perturb(z,w,eta,beta,seed)
use paran
use commk
use lyap
use commphy
use cudafor
implicit none

complex(8), dimension(NX2P1,NY) :: w,z
real(8) k2,k,fac,eta,beta,scax,scay,scra,scab
real(8) i0,j0,k0,alx,aly,s2,s1
integer :: i,j,seed,uu
real :: rann

i0 = rann(seed)*xlx
j0 = rann(seed)*xly

! s1=xlx/2.d0
! if (xlx.gt.xly) then
! 	s1=xly/2.d0
! endif
! 
! do uu=1,100
! 	k0 = rann(seed)*xlx
! 	alx= mod( i0 + xlx/2.d0 * dcos(k0),xlx)
! 	aly= mod( j0 + xly/2.d0 * dsin(k0),xly)
! 	if (alx.lt.0) then
! 		alx=alx+xlx
! 	endif
! 	if (aly.lt.0) then
! 		aly=aly+xly
! 	endif
! 	
! 	s2 = dsqrt((alx-i0)**2.d0+(aly-j0)**2.d0)
! 	if ((s2.gt.(0.99*s1)).and.(s2.lt.(1.01*s1))) exit
! enddo

! write(6,*),"i0x,j0y, = ",real(i0),real(j0)
! call flush(6)
! write(6,*),"alx,aly, = ",real(alx),real(aly),real(s2)
! call flush(6)
! write(6,*),"s2 = ",real(s2)
! call flush(6)

!$acc kernels present(eed,w)
eed = 0.d0
w = 0.d0
!$acc end kernels

!$acc parallel loop collapse(2) present(w,fkxs,fkys,fkx,fky,eed)
do 10 j=1,NY
do 10 i=1,NX2P1
  if ((i.eq.1).and.(j.eq.1)) goto 10
  fac=1.d0
  if ((i.eq.1).or.(i.eq.NX2P1)) fac=0.5d0
  k2=fkxs(i) + fkys(j)
  k = sqrt(k2)
  scax = fkx(i)*i0  + fky(j)*j0
!   scay = fkx(i)*alx + fky(j)*aly
  scab = dexp(-0.5d0*(k*eta)**2)*k**(1-beta)
!   w(i,j) = (cdexp(-csqrt(-1.d0)*scax) - cdexp(-csqrt(-1.d0)*scay) ) * scab
    w(i,j) = (cdexp(-csqrt(-1.d0)*scax)) * scab
!   w(i,j) = scab
  scra=0.5d0*cdabs( w(i,j) )**2 
  !$acc atomic
  eed=eed+fac*scra      ! Z[dw](t)
10 continue
!$acc end parallel loop

! Norm delta0^2/2
!$acc parallel loop collapse(2) present(w,z,eed)
do j=1,NY
	do i=1,NX2P1
		w(i,j) = z(i,j) + delta0 * w(i,j) / dsqrt(2.d0*eed)
	enddo
enddo
!$acc end parallel loop

return
end
