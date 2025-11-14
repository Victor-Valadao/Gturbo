! Computation of important physical operators, including:
! nonlinear term         -> nlt()
! velocity field         -> vel()
! physical space product -> produ()
! circular truncation f  -> trunc()

!=================================================

subroutine vel(z,u,v)
use paran
use cufft
use openacc
use cudafor
use stream
use plan
use commk
use commphy
implicit none

real(8), dimension(2,NX2P1,NY) :: z,u,v
real(8) k2,ka,kx,ky
integer i,j

!$acc kernels present(u,v)
u=0.d0
v=0.d0
!$acc end kernels

!$acc parallel loop collapse(2) present(u,v,z,fkxs,fkys,fkx,fky)
do 10 j=1,NY
do 10 i=1,NX2P1
if ((i.eq.1).and.(j.eq.1)) goto 10
	k2=fkxs(i)+fkys(j)
	ka=k2**(alpv/2.d0)
	kx=fkx(i)/ka
	ky=fky(j)/ka
	u(1,i,j)=-ky*z(2,i,j)
	u(2,i,j)=+ky*z(1,i,j)
	v(1,i,j)=+kx*z(2,i,j)
	v(2,i,j)=-kx*z(1,i,j)
10 continue
!$acc end parallel

return
end

!=================================================

subroutine nlt(z,u,v)
use paran
use cufft
use openacc
use cudafor
use stream
use plan
use commk
use commphy
implicit none

real(8), intent(inout ), dimension(2,NX2P1,NY) :: z
real(8), dimension(2,NX2P1,NY) :: u, v
real(8) k2,ka,kx,ky
integer i,j

! Calculate velocity field

call vel(z,u,v)

! IFFT to Physical space u,v,z. Some improvement can be done here
! to perform all FFT's at the same time since they are independent, 
! but it will consume more memory
call fft_inv(z, plan_inv)
call fft_inv(v, plan_inv)
call fft_inv(u, plan_inv)

! Calculate the product \omega v_i in the real space and store in v_i
call produ(z,u,v)

! FFT back to Fourier space. Same paralel improvement can be
! done here at the cost of memory
call fft_dir(u,plan_dir)
call fft_dir(v,plan_dir)

! Calculate Nonlinear term in Fourier space div(\omega v_i)
!$acc parallel loop collapse(2) present(z, u, v, fkx, fky)
do j=1,NY
do i=1,NX2P1
	z(1,i,j)=+fkx(i)*u(2,i,j)+fky(j)*v(2,i,j)
	z(2,i,j)=-fkx(i)*u(1,i,j)-fky(j)*v(1,i,j)
end do
end do
!$acc end parallel loop

if (truncate.eq.1) call trunc(z) ! dealiasing

return
end

!=================================================

subroutine produ(z,u,v)
use paran
use cudafor
implicit none

real(8), intent(in    ), dimension(NXP2,NY) :: z
real(8), intent(inout   ), dimension(NXP2,NY) :: u
real(8), intent(inout   ), dimension(NXP2,NY) :: v
integer :: i,j

!$acc parallel loop collapse(2) present(u, v, z) 
do j=1,NY
do i=1,NXP2
  if(i.le.NX) then
    u(i,j)=z(i,j)*u(i,j)
    v(i,j)=z(i,j)*v(i,j)
  else
    u(i,j)=0.d0
    v(i,j)=0.d0
  end if
end do
end do
!$acc end parallel

return
end

!=================================================

subroutine trunc(z)
use paran
use commk
implicit none

real(8), dimension(2,NX2P1,NY) :: z
real(8) :: k1,k2
integer :: i,j

!$acc parallel loop collapse(2) present(z, fkxs, fkys)
do j=1,NY
 do i=1,NX2P1
  k1=fkys(j)/fkysmax
  k2=k1+fkxs(i)/fkxsmax
  if (k2.gt.rtrunc) then
   z(1,i,j)=0.d0
   z(2,i,j)=0.d0
  end if
 end do
end do
!$acc end parallel

! Subtract a possible existing mean z
!$acc kernels present(z)
z(1,1,1)=0.d0
z(2,1,1)=0.d0
!$acc end kernels

return
end