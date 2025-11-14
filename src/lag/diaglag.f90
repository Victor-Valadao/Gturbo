! Write Lagrangian observables, including:
! writeftle()  <- writes './files/ftle*.*'
! writelyap()  <- writes './files/lyapunov.*'
! writelag()   <- writes './fields/Xp.*' and './fields/Dp.*'

!=================================================

subroutine writeftle
use modlag
use commphy
use commsim
implicit none
integer ifr,i

do i=1,NP
write(22,*)real(ftle(i,1))
write(23,*)real(ftle(i,2))
enddo

call flush(22)
call flush(23)

return
end

!=================================================

subroutine writelyap(t0)
use modlag
implicit none
integer :: iv
real(8) :: t0

write(17,99)real(t0),(real(lyap(iv)),iv=1,NV)
99 format(3g)
call flush(17)

return 
end

!=================================================

subroutine writelag(xp,dxp)
use paran
use modlag
implicit none

real(8) :: xp(2,NP),dxp(NV,2,NP)

write(40,*)xp
call flush(40)

write(41,*)dxp
call flush(41)  
    
return
end

!=================================================

subroutine write_int_vel_pos(z,u,v,dux,duy,dvx,xp,dxp,t0)
use paran
use cufft
use openacc
use cudafor
use stream
use plan
use commk
use commphy
use modlag
implicit none

real(8), intent(inout ), dimension(2,NX2P1,NY) :: z
real(8), dimension(2,NX2P1,NY) :: u, v, dux, duy, dvx
real(8), dimension(NV,2,NP) :: dxp
real(8), dimension(2,NP) :: xp,p1,p2,v1,v2
real(8) :: t0,rebox2
integer i,ip

call vel(z,u,v)
call grads(u,v,dux,duy,dvx)

call fft_inv(z, plan_inv)
call fft_inv(v, plan_inv)
call fft_inv(u, plan_inv)
call fft_inv(dux, plan_inv)
call fft_inv(duy, plan_inv)
call fft_inv(dvx, plan_inv)

p2=0.d0
p1=0.d0
v2=0.d0
v1=0.d0
!$acc enter data copyin(p1,p2,v1,v2)

!it would be good to save the velocity gradients in the particle position

!$acc parallel loop collapse(2) present(xp,p1,p2,v1,v2,dxp)
do ip=1,NP
do i=1,nv
	p1(i,ip)=xp(i,ip)
	p2(i,ip)=xp(i,ip)+dxp(1,i,ip)
	
	v1(i,ip)=xp(i,ip)
	v2(i,ip)=xp(i,ip)+dxp(1,i,ip)
enddo
enddo
!$acc end parallel
 
!$acc parallel loop collapse(1) present(p2,v2)
do ip=1,NP
 p2(1,ip)=rebox2(p2(1,ip),xlx)
 p2(2,ip)=rebox2(p2(2,ip),xly)
 v2(1,ip)=rebox2(v2(1,ip),xlx) ! this does not exist right?
 v2(2,ip)=rebox2(v2(2,ip),xly)
end do
!$acc end parallel

call interpolate(u,v,dux,duy,dvx,v1,dxp)
call interpolate(u,v,dux,duy,dvx,v2,dxp)

!$acc update host(p1,p2,v1,v2)

do i=1,NP
write(19,99)real(t0),real(p1(1,i)),real(p1(2,i)),real(v1(1,i)),real(v1(2,i)),real(p2(1,i)),real(p2(2,i)),real(v2(1,i)),real(v2(2,i))
99 format(9g)
enddo

call flush(19)
!$acc exit data delete(p1,p2,v1,v2) finalize
return 
end

! --------------------------------------------------------------
! Reboxing function
function rebox2(x,xl)
implicit none
real(8) :: x,xl,rebox2

rebox2=x
if (rebox2.lt.0.d0) then
 rebox2=rebox2+xl
end if

if (rebox2.ge.xl) then
 rebox2=rebox2-xl
end if

return
end 

