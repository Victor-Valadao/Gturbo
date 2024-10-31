subroutine trunc(z)
use paran
use commk
implicit none

complex*16 z(NX2P1,NY)
real(8) :: k1,k2
integer :: i,j

!$acc parallel loop collapse(2) present(z, fkxs, fkys)
do j=1,NY
 do i=1,NX2P1
  k1=fkys(j)/fkysmax
  k2=k1+fkxs(i)/fkxsmax
  if (k2.gt.rtrunc) then
   z(i,j)=(0.d0,0.d0)
  end if
 end do
end do
!$acc end parallel

! Substract a possible existing mean z
!$acc kernels present(z)
z(1,1)=(0.d0,0.d0)
!$acc end kernels

return
end
! --------------------------------------------------------------	
! Sub to print memory

subroutine show_mem(fre, tot, step)
use cudafor
implicit none

integer :: ierr,step
integer(kind=cuda_count_kind) :: fre,tot
real :: r

ierr=cudaMemGetInfo( fre, tot )
r=(tot-fre)/2.0**30

write(6,"(A25,I4,A3,F12.9,A4,F6.3)")"Alloc mem (Gb) in step", step, " = ",r," of ",tot/2.0**30
call flush(6)

tot=0
fre=0

end