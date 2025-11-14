subroutine strucs(z,u,v)
use paran
use cufft
use openacc
use cudafor
use stream
use plan
use commsim
use commphy
implicit none

real(8), intent(in ), dimension(NXP2,NY) :: z,u,v
real(8), dimension(nscale,8)   :: szz,sul,sut,slz,stz,s12
real(8) :: dzz, dut, dul, stat
integer :: is, i, j, ll, ii, jj

real(8), dimension(6,8) :: temp  ! Renamed from t to temp
temp = 0.d0
stat = 0.d0

!$acc enter data copyin(temp,stat,szz,sul,sut,slz,s12)

!$acc kernels present(szz,sul,sut,slz,s12)
szz = 0.d0
sul = 0.d0
sut = 0.d0
slz = 0.d0
s12 = 0.d0
!$acc end kernels

do is = 1, nscale

 !$acc kernels present(temp,stat)
  temp = 0.d0
  stat = 0.d0
 !$acc end kernels

 !$acc parallel loop collapse(2) present(z,u,v,temp,stat)
  do j = 1, NY, nskip
    do i = 1, NX, nskip

      ii = i + iss(is); if (ii > NX) ii = ii - NX
      dzz = z(ii,j) - z(i,j)
      dul = u(ii,j) - u(i,j)
      dut = v(ii,j) - v(i,j)
      do ll = 1,8
        !$acc atomic
        temp(1,ll) = temp(1,ll) + abs(dzz)**(ll)
        !$acc atomic
        temp(2,ll) = temp(2,ll) + abs(dul)**(ll)
        !$acc atomic
		temp(3,ll) = temp(3,ll) + abs(dut)**(ll)
		!$acc atomic
        temp(4,ll) = temp(4,ll) + dul*abs(dzz)**(ll-1)
!         temp(5,ll) = temp(5,ll) + dut*(dzz)**ll
		!$acc atomic
		temp(6,ll) = temp(6,ll) + abs( dul * dzz**(2) )**(dble(ll)/3.d0)

!         if (mod(ll-1,3) .eq. 0) then
!           temp(6,ll) = temp(6,ll) + (dul*dzz*dzz)**(dble(ll-1)/3.d0)
!         else
!           temp(6,ll) = temp(6,ll) + abs(dul*dzz*dzz)**(dble(ll-1)/3.d0)
!         end if
      end do
      !$acc atomic
      stat = stat + 1.d0

      jj = j + iss(is); if (jj > NY) jj = jj - NY
      dzz = z(i,jj) - z(i,j)
      dul = v(i,jj) - v(i,j)
      dut = u(i,jj) - u(i,j)
      do ll = 1,8
		!$acc atomic
        temp(1,ll) = temp(1,ll) + abs(dzz)**(ll)
        !$acc atomic
        temp(2,ll) = temp(2,ll) + abs(dul)**(ll)
        !$acc atomic
		temp(3,ll) = temp(3,ll) + abs(dut)**(ll)
		!$acc atomic
        temp(4,ll) = temp(4,ll) + dul*abs(dzz)**(ll-1)
!         temp(5,ll) = temp(5,ll) + dut*abs(dzz)**(ll-1)
		!$acc atomic
		temp(6,ll) = temp(6,ll) + abs( dul * dzz**(2) )**(dble(ll)/3.d0)
		
!         if (mod(ll-1,3) .eq. 0) then
!           temp(6,ll) = temp(6,ll) + (dul*dzz*dzz)**(dble(ll-1)/3.d0)
!         else
!           temp(6,ll) = temp(6,ll) + abs(dul*dzz*dzz)**(dble(ll-1)/3.d0)
!         end if
      end do
      !$acc atomic
      stat = stat + 1.d0

    end do
  end do
 !$acc end parallel loop

!  !$acc update host(temp,stat)

  ! Normalize and store
  !$acc parallel loop present(temp,stat,szz,sul,sut,slz,s12)
  do ll = 1, 8
    szz(is,ll) = temp(1,ll) / stat
    sul(is,ll) = temp(2,ll) / stat
    sut(is,ll) = temp(3,ll) / stat
    slz(is,ll) = temp(4,ll) / stat
!     stz(is,ll) = temp(5,ll) / stat
    s12(is,ll) = temp(6,ll) / stat
  end do
end do

!$acc update host(szz,sul,sut,slz,s12,stat)
! Write results
do j = 1, nscale
  write(40,99)real(stat), real(iss(j)), &
       (szz(j,ll), ll=1,8), &
       (sul(j,ll), ll=1,8), &
       (sut(j,ll), ll=1,8), &
       (slz(j,ll), ll=1,8), &
       (s12(j,ll), ll=1,8)
  99 format(42e14.6)
end do

! do j = 1, nscale
!   write(40,99)real(stat), real(iss(j)), &
!        (szz(j,ll), ll=1,10), &
!        (sul(j,ll), ll=1,10), &
!        (sut(j,ll), ll=1,10), &
!        (slz(j,ll), ll=1,10), &
!        (stz(j,ll), ll=1,10), &
!        (s12(j,ll), ll=1,10)
!   99 format(62e14.6)
! end do

!$acc exit data delete(temp,stat,szz,sul,sut,slz,s12)

return
end
