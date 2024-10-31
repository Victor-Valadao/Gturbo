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