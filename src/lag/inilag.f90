! Lagrangian initialization
! load / initialize positions and perturbation vectors

! Initialization of the Lagrangian part
subroutine inilag(xp,dxp,ifr,seed)
use paran
use commphy
use commsim
use modlag
use commforc
implicit none

real(8), intent(inout) :: xp(2,NP), dxp(NV,2,NP)
! real(4) :: Ap(2,NP), Axp(NV,2,NP)
real(8) :: pi2,time0,prod,a
parameter(pi2=2.d0*3.14159265358979d0)
integer i,j,ip,im,ifr,seed
character*16 nome
real :: rann

!Fix interpolation scale
scalax=dble(NX)/xlx
scalay=dble(NY)/xly
ftle=0.d0
lyap=0.d0
time0=0.d0
co=0.d0

! Lyapunov initialization
if (ifr==0) then
 do ip=1,NP
 
 ! particles positions x=[0,2pi)
 xp(1,ip)=pi2*rann(seed)
 xp(2,ip)=pi2*rann(seed)
 
  ! random orthonormal vectors
  a=2.d0*(rann(seed)-0.5)
  dxp(1,1,ip)=a
  dxp(1,2,ip)=sqrt(1-a**2)
  dxp(2,1,ip)=sqrt(1-a**2)
  dxp(2,2,ip)=-a
 !  ! random vector
!   do i=1,NV
!    dxp(i,1,ip)=rann(seed)-0.5
!    dxp(i,2,ip)=rann(seed)-0.5
!   enddo
!   
!   ! Gram-Schimdt to orthogonalize it
!   do i=1,nv  
!    do j=1,i-1
!     prod = dxp(i,1,ip)*dxp(j,1,ip) + dxp(i,2,ip)*dxp(j,2,ip)
!     dxp(i,1,ip) = dxp(i,1,ip) - prod*dxp(j,1,ip)
!     dxp(i,2,ip) = dxp(i,2,ip) - prod*dxp(j,2,ip)
!    end do
!    
!    ! Normalization 
!    prod=dsqrt(dxp(i,1,ip)*dxp(i,1,ip)+dxp(i,2,ip)*dxp(i,2,ip))
!    dxp(i,1,ip)=dxp(i,1,ip)/prod
!    dxp(i,2,ip)=dxp(i,2,ip)/prod
!   end do
 end do
  write(nome,"('./fields/Dp.',i3.3)")ifr
  open(unit=70,file=nome,form='unformatted')
  write(70)dxp
  close(70)
  
  write(nome,"('./fields/Xp.',i3.3)")ifr
  open(unit=70,file=nome,form='unformatted')
  write(70)xp
  close(70)
 
else ! Load previous simulation
  write(nome,"('./fields/Dp.',i3.3)")ifr
  open(unit=70,file=nome,form='unformatted')
  read(70)dxp
  close(70)
  
  write(nome,"('./fields/Xp.',i3.3)")ifr
  open(unit=70,file=nome,form='unformatted')
  read(70)xp
  close(70)
endif   

! xp=dble(Ap)
! dxp=dble(Axp)  

write(6,*)' '
write(6,*)' ---------------------------------------------------'
if (oitrp.eq.1) then
  !1st order
  write(6,*)' Lagrangian Lyapunov with 1st ord interpolation'
  co(1,1,1)= 1.d0
  co(1,1,2)=-1.d0
  co(1,2,2)= 1.d0
else if (oitrp.eq.2) then
  write(6,*)' Lagrangian Lyapunov with 2st ord interpolation'
  !2nd order
  co(2,1,1)= 1.0d0
  co(2,1,2)=-1.5d0
  co(2,1,3)= 0.5d0

  co(2,2,2)= 2.d0
  co(2,2,3)=-1.d0

  co(2,3,2)=-0.5d0
  co(2,3,3)= 0.5d0
else
  write(6,*)' Lagrangian order invalid, falling back to 1st'
  oitrp=1
  co(1,1,1)= 1.d0
  co(1,1,2)=-1.d0
  co(1,2,2)= 1.d0
end if
write(6,*)' Particle number = ',NP
! write(6,*)' OITRP = ',oitrp
write(6,*)' ---------------------------------------------------'
write(6,*)' '

call flush(6)

return
end
