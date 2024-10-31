subroutine forcing(z,seed)
!$acc routine(rann) seq
use paran
use commforc
use commk
use commphy
implicit none

real(8) :: pi2
parameter(pi2=2.d0*3.14159265358979d0)
real(8) :: te,cc,cs
real(8), intent(inout ), dimension(2,NX2P1,NY) :: z
integer :: i,j,k,seed
real :: rann

!$acc parallel loop seq present(z,jc,ff,ikxf,ikyf,seed,sdt)
do k=1,nf
 i=ikxf(k)
 j=ikyf(k)
 
 te=pi2*rann(seed)
 cc=cos(te)
 cs=sin(te)

 z(1,i,j)=z(1,i,j)+sdt*famp*ff(1,k)*cc
 z(2,i,j)=z(2,i,j)+sdt*famp*ff(2,k)*cs

 if ((i.eq.1).or.(i.eq.NX2P1)) then
	z(1,i,jc(j))=z(1,i,jc(j))+sdt*famp*ff(1,k)*cc
	z(2,i,jc(j))=z(2,i,jc(j))-sdt*famp*ff(2,k)*cs
 end if
end do
!$acc end parallel loop

return
end

! Subroutine to generate the random seed for the forcing

real function rann(irand)
!$acc routine seq  
implicit none
integer irand
integer mask

mask = '7FFFFFFF'X
irand = iand((69069*irand + 1),mask)
rann = dble(irand) / dble(mask)

return
end