! Subroutines to read field in real/spectral space with single precision

subroutine readf(z,ifr)
  use paran
  integer i, j,ifr
  real*8 z(NXP2,NY)
real*4 w(NX,NY)
character*40 nome

write(nome,"('./fields/w.',i3.3)")ifr
open(unit=90,file=nome,form='unformatted',action='read')
 read(90)w
close(90)

do j=1,NY
do i=1,NX
 z(i,j)=dble(w(i,j))
end do
end do

! Fill last lines with 0 to perform the FFT correctly
do j=1,NY
do i=NX+1,NXP2
 z(i,j)=0.d0
end do
end do

return
end

subroutine readfk(u,ifr)
use paran
implicit none

integer :: ifr,i,j
real(8), dimension(NXP2,NY) :: u
character*40 nome

write(nome,"('./Diag/fields/wk.',i8.8)")ifr
open(unit=90,file=nome,form='unformatted')

read(90)u
close(90)

return
end