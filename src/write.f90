! Subroutines to write field in real/spectral space with single precision

subroutine writef(u,ifr)
use paran
use commphy
implicit none

integer :: ifr,i,j
real(8), dimension(NXP2,NY) :: u
real(4), dimension(NX,NY) :: w
character*40 nome

write(nome,"('./fields/w.',i3.3)")ifr
open(unit=90,file=nome,form='unformatted')

do j=1,NY
do i=1,NX
 w(i,j)=real(u(i,j))
end do
end do
write(90)w
close(90)

return
end

!___________________________

subroutine writef2(u,ifr,jfr)
use paran
use commphy
implicit none

integer :: ifr,jfr,i,j
real(8), dimension(NXP2,NY) :: u
real(4), dimension(NX,NY) :: w
character*40 nome

write(nome,"('./fields/w.',i3.3,'.',i3.3)")ifr,jfr
open(unit=90,file=nome,form='unformatted')

do j=1,NY
do i=1,NX
 w(i,j)=real(u(i,j))
end do
end do
write(90)w
close(90)

return
end

subroutine writefk(u,ifr)
use paran
implicit none

integer :: ifr,i,j
real(8), dimension(NXP2,NY) :: u,v
character*40 nome

write(nome,"('./Diag/fields/wk.',i8.8)")ifr
open(unit=90,file=nome,form='unformatted')


do j=1,NY
do i=1,NX
 v(i,j)=real(u(i,j))
end do
end do
write(90)v
close(90)

return
end