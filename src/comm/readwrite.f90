! Read / Write subroutines on single and double precision, including:
! readf()   <- reads  './fields/w.*'
! writef()  <- writes './fields/w.*'
! writef2() <- writes './fields/w.*.*' <- intermediate saves
! readp()   <- reads  './fields/p.*'
! writep()  <- writes './fields/p.*'
! writep2() <- writes './fields/p.*.*' <- intermediate saves

! Also includes subroutine to report 
! memory consumption inside the GPU

!=================================================

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

!=================================================

subroutine readf(z, ifr)
use paran
use commphy
implicit none

integer :: i, j, ifr
real(8), dimension(NXP2, NY) :: z
real(4), dimension(NX, NY) :: ws
real(8), dimension(NX, NY) :: wd
character(len=40) :: nome

write(nome, "('./fields/w.', i3.3)") ifr
open(unit=90, file=nome, form='unformatted', action='read')

if (prec == 4) then
 read(90) ws
 do j=1,NY
 do i=1,NX
  z(i,j)=dble(ws(i,j))
 end do
 end do

else if (prec == 8) then
 read(90) wd
 do j=1,NY
 do i=1,NX
  z(i,j) = wd(i,j)
 end do
 end do

else
 write(6,*)'Error: unsupported precision (prec=', prec, ')'
 call flush(6)
end if

close(90)

! Fill last lines with 0 to perform the FFT correctly
do j=1,NY
do i=NX+1,NXP2
 z(i,j)=0.d0
end do
end do

return
end

!=================================================

subroutine writef(u, ifr)
use paran
use commphy
implicit none

integer :: ifr, i, j
real(8), dimension(NXP2, NY) :: u
real(4), dimension(NX, NY) :: ws
real(8), dimension(NX, NY) :: wd
character(len=40) :: nome

write(nome, "('./fields/w.', i3.3)") ifr
open(unit=90, file=nome, form='unformatted')

if (prec == 4) then
 do j=1,NY
 do i=1,NX
  ws(i,j) = real(u(i,j), kind=4)
 end do
 end do
 write(90) ws

else if (prec == 8) then
 do j=1,NY
 do i=1,NX
  wd(i,j) = u(i,j)
 end do
 end do
 write(90) wd

else
 write(6,*)'Error: unsupported precision (prec=', prec, ')'
 call flush(6)
end if

close(90)
return
end

!=================================================

subroutine writef2(u,ifr,jfr)
use paran
use commphy
implicit none

integer :: ifr,jfr,i,j
real(8), dimension(NXP2, NY) :: u
real(4), dimension(NX, NY) :: ws
real(8), dimension(NX, NY) :: wd
character(len=40) :: nome

write(nome,"('./fields/w.',i3.3,'.',i3.3)")ifr,jfr
open(unit=90,file=nome,form='unformatted')

if (prec == 4) then
 do j=1,NY
 do i=1,NX
  ws(i,j) = real(u(i,j), kind=4)
 end do
 end do
 write(90) ws

else if (prec == 8) then
 do j=1,NY
 do i=1,NX
  wd(i,j) = u(i,j)
 end do
 end do
 write(90) wd

else
 write(6,*)'Error: unsupported precision (prec=', prec, ')'
 call flush(6)
end if

close(90)
return
end

!=================================================

subroutine readp(z, ifr)
use paran
use commphy
implicit none

integer :: i, j, ifr
real(8), dimension(NXP2, NY) :: z
real(4), dimension(NX, NY) :: ws
real(8), dimension(NX, NY) :: wd
character(len=40) :: nome

write(nome, "('./fields/p.', i3.3)") ifr
open(unit=90, file=nome, form='unformatted', action='read')

if (prec == 4) then
 read(90) ws
 do j=1,NY
 do i=1,NX
  z(i,j)=dble(ws(i,j))
 end do
 end do

else if (prec == 8) then
 read(90) wd
 do j=1,NY
 do i=1,NX
  z(i,j) = wd(i,j)
 end do
 end do

else
 write(6,*)'Error: unsupported precision (prec=', prec, ')'
 call flush(6)
end if

close(90)

! Fill last lines with 0 to perform the FFT correctly
do j=1,NY
do i=NX+1,NXP2
 z(i,j)=0.d0
end do
end do

return
end

!=================================================

subroutine writep(u, ifr)
use paran
use commphy
implicit none

integer :: ifr, i, j
real(8), dimension(NXP2, NY) :: u
real(4), dimension(NX, NY) :: ws
real(8), dimension(NX, NY) :: wd
character(len=40) :: nome

write(nome, "('./fields/p.', i3.3)") ifr
open(unit=90, file=nome, form='unformatted')

if (prec == 4) then
 do j=1,NY
 do i=1,NX
  ws(i,j) = real(u(i,j), kind=4)
 end do
 end do
 write(90) ws

else if (prec == 8) then
 do j=1,NY
 do i=1,NX
  wd(i,j) = u(i,j)
 end do
 end do
 write(90) wd

else
 write(6,*)'Error: unsupported precision (prec=', prec, ')'
 call flush(6)
end if

close(90)
return
end

!=================================================

subroutine writep2(u,ifr,jfr)
use paran
use commphy
implicit none

integer :: ifr,jfr,i,j
real(8), dimension(NXP2, NY) :: u
real(4), dimension(NX, NY) :: ws
real(8), dimension(NX, NY) :: wd
character(len=40) :: nome

write(nome,"('./fields/p.',i3.3,'.',i3.3)")ifr,jfr
open(unit=90,file=nome,form='unformatted')

if (prec == 4) then
 do j=1,NY
 do i=1,NX
  ws(i,j) = real(u(i,j), kind=4)
 end do
 end do
 write(90) ws

else if (prec == 8) then
 do j=1,NY
 do i=1,NX
  wd(i,j) = u(i,j)
 end do
 end do
 write(90) wd

else
 write(6,*)'Error: unsupported precision (prec=', prec, ')'
 call flush(6)
end if

close(90)
return
end

