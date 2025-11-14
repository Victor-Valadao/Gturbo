! Defines wave vectors, linear integrator operator and others.

subroutine defk
use paran
use commphy
use commk
use commforc
use pass
implicit none

real(8) :: fk2
real(8) :: scra,km2
integer :: i,j

! Conjugate variables: jc
jc(1)=1
do j=2,NY
 jc(j)=NY+2-j
end do

! Wavenumber increments
dkx=2.0d0*3.14159265358979323846d0/xlx
dky=2.0d0*3.14159265358979323846d0/xly

! Wavenumbers
do i=1,NX2P1
 fkx(i)=dble(i-1)*dkx
end do
do j=1,NY2P1
 fky(j)=dble(j-1)*dky
end do
do j=NY2P2,NY
 fky(j)=dble(j-NY-1)*dky
end do

! Wavenumbers squared 
do i=1,NX2P1
 fkxs(i)=fkx(i)*fkx(i)
end do
do j=1,NY
 fkys(j)=fky(j)*fky(j)
end do

! Dissipative terms are exponentially integrated e**D*dt/2
! D=nu*k2**alpnu+mu*k2**alpmu

emk(1,1) = 0.d0
epk(1,1) = 0.d0
do 15 j=1,NY
do 15 i=1,NX2P1
 if ((i.eq.1).and.(j.eq.1)) goto 15
 fk2=fkxs(i)+fkys(j)
 scra=nu*fk2**alpnu+mu*fk2**alpmu
 emk(i,j)=dexp(-0.5d0*dt*scra)  
 scra=nup*fk2**alpnup+mup*fk2**alpmup
 epk(i,j)=dexp(-0.5d0*dt*scra)
15 continue

! Maximum wavenumbers
fkxsmax=fkxs(NX2P1)
fkysmax=fkys(NY2P1)
rtrunc=rtrunc*rtrunc

return
end