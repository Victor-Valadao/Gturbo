use cudafor
use openacc
use cufft
use iso_c_binding , only: C_PTR
use paran
use stream
use plan
use commsim
use commphy
use commk
use commforc
use modlag
use time_keeper
implicit none

real(8), dimension(NS) :: z,znl,zne,u,v,dux,duy,dvx
real(8), dimension(NV,2,NP) :: dp,dpnl,dpne
real(8), dimension(2,NP) :: xp,xpnl,xpne
integer :: ipas,ifr,jfr,istat,seed,seedp,ierr,err_dir,err_inv,i,j
integer :: startTime, stopTime, trate
integer(kind=cuda_count_kind) :: fre,tota
real(8) :: s,t0,et_s
character*40 nome

call startclock()
! call system_clock(startTime)

! Shows the amount of memory consumed in the GPU
call show_mem(fre,tota,0)

! InvFFT Complex to Real double precision
err_inv = err_inv + cufftPlan2D(plan_inv,NY,NY,CUFFT_Z2D)
istat = cudaStreamCreate(stream1)
err_inv = err_inv + cufftSetStream(plan_inv,stream1)

! DirFFT Real to Complex double precision
err_dir = err_dir + cufftPlan2D(plan_dir,NY,NY,CUFFT_D2Z)
istat = cudaStreamCreate(stream2)
err_dir = err_dir + cufftSetStream(plan_dir,stream2)

! Reads if the simulation is a continuation of a previous one
open(unit=1,file='./curframe.dat')
 read(1,*)ifr,t
close(1)

! Reads the random seed
write(nome,"('./files/seed.',i3.3)")ifr
open(unit=2,file=nome)
 read(2,*)seed
close(2)

seedp=seed
! Initializes wave number, forcing, initial field z (real space)
call inieuler(z, ifr, seedp)

! Opening files to save global quantities and useful information
write(nome,"('./files/lyapunov.',i3.3)")ifr
Open(unit=17,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/global.',i3.3)")ifr
Open(unit=18,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/local.',i3.3)")ifr
Open(unit=19,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/spectra.',i3.3)")ifr
Open(unit=20,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/ftle1.',i3.3)")ifr
open(unit=22,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/ftle2.',i3.3)")ifr
open(unit=23,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/fluxes.',i3.3)")ifr
Open(unit=30,file=nome,form='formatted',status='unknown',access='append')

seedp=seed
call inilag(xp,dp,ifr,seedp)

! Shows the amount of memory consumed in the GPU
call show_mem(fre,tota,1)

! ========================================================

! Copying data to the GPU data

!$acc enter data copyin(emk)                  !<-- calc in defk called in inieuler
!$acc enter data copyin(ikxf,ikyf,ff)         !<-- calc in iniforcing called in inieuler
!$acc enter data copyin(fkx,fky,fkxs,fkys,jc) !<-- calc in defk called in inieuler
!$acc enter data copyin(dt,dt2,dt3,dt6,sdt)   !<-- calc in inieuler/loaded
!$acc enter data copyin(seed)                 !<-- read in main used in forcing/rann
!$acc enter data copyin(znl,zne,u,v)          !<-- allocate memory for znl, u, v
!$acc enter data copyin(dux,duy,dvx)
!$acc enter data copyin(co,lyap,ftle)
!$acc enter data copyin(xp,xpnl,xpne)
!$acc enter data copyin(dp,dpnl,dpne)

!$acc enter data copyin(z)           !<-- pass the initial field to the GPU

! Shows the amount of memory consumed in the GPU
call show_mem(fre,tota,2)

! DirFFT the initial field + dealiasing
call fft_dir(z,plan_dir)
if (truncate.eq.1) then
 call trunc(z)
end if

jfr=1
! =============================================================
!  Time evolution <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
do ipas=1,npas
  t=t+dt
  
  if (rko.eq.4) then
  	call step4(z,znl,zne,u,v,dux,duy,dvx,xp,dp,xpnl,dpnl,xpne,dpne)
  else
  	call step(z,znl,u,v,dux,duy,dvx,xp,dp,xpnl,dpnl)
  end if
  
  if (mod(ipas,imix).eq.0) then  ! Write diagnostics each imix
    !$acc kernels present(z,znl,xp,xpnl,dp,dpnl)
	znl = z
	xpnl= xp
	dpnl= dp
	!$acc end kernels
	call nlt(znl,u,v,dux,duy,dvx,xpnl,dpnl)

  	call spectrum(z,ipas)        ! Write Energy and Enstrophy spectra
    call fluxes(z,znl,ipas)      ! Write Energy and Enstrophy fluxes
    call glob(z,ipas)            ! Write Dissipative terms
  endif
  
  if ((mod(ipas,iout).eq.0).and.(ipas.ne.npas)) then  ! Save vorticity each iout
    !$acc kernels present(z,znl) ! Copy z to znl
	znl = z
	!$acc end kernels
	
	call fft_inv(znl , plan_inv)  ! FFt to save in the physical space
	!$acc update host(znl,xp,dp) 
  	
  	write(nome,"('./fields/Dp.',i3.3,'.',i3.3)")ifr,jfr
	open(unit=70,file=nome,form='unformatted')
	write(70)dp
	close(70)

	write(nome,"('./fields/Xp.',i3.3,'.',i3.3)")ifr,jfr
	open(unit=70,file=nome,form='unformatted')
	write(70)xp
	close(70)
  	call writef2(znl ,ifr,jfr)        ! Save field
  	jfr=jfr+1
  	
  	write(6,*)"Saved w field, seed and curframe at timestep = ",ipas
  	call flush(6)
  	
  endif
  
	! ----------------------------------------------------------
	! Lyapunov
	if (mod(ipas,nlyal).eq.0) then 
	 ! ! 
! 	 !$acc kernels present(z,znl,xp,xpnl,dp,dpnl)
! 	 znl = z
! 	 xpnl= xp
! 	 dpnl= dp
! 	 !$acc end kernels
!      call write_int_vel_pos(znl,u,v,dux,duy,dvx,xpnl,dpnl,t)
	
	
 	 call rescale(dp)
	 !$acc update host(lyap,ftle)

! 	 !$acc update host(dp)
! 	 call rescale(dp)
! 	 !$acc update device(dp)
	 
	 t0=dble(dt*nlyal)
	 lyap=lyap/t0
	 ftle=ftle/t0
	 call writelyap(t)
	 call writeftle
	 write(6,*)"Nt = ",ipas,", ly = ",real(lyap(1))," and ",real(lyap(2))
! 	 call show_mem(fre,tota,2)
	call flush(6)
	end if
	! ----------------------------------------------------------
  
  ! ! Add forcing contribution
  call forcing(z,seed)
  
end do
! --------------------------------------

ifr=ifr+1
! Final write
!$acc update host(z)
call fft_inv(z, plan_inv)

!$acc update host(seed)
write(nome,"('./files/seed.',i3.3)")ifr
open(unit=2,file=nome)
write(2,*)seed
close(2)

! Export data to curframe.dat to continue, if you want
open(unit=1,file='./curframe.dat')
 write(1,*)ifr,t
close(1)

close(17)
close(18)
close(19)
close(20)
close(22)
close(23)
close(30)

! ========================================================

! Exit GPU data Freeing the memory

!$acc exit data delete(emk) finalize
!$acc exit data delete(ikxf,ikyf,ff) finalize
!$acc exit data delete(fkx,fky,fkxs,fkys,jc) finalize
!$acc exit data delete(dt,dt2,dt3,dt6,sdt) finalize
!$acc exit data delete(seed) finalize
!$acc exit data delete(co,lyap,ftle) finalize
!$acc exit data delete(u,v,znl,zne) finalize
!$acc exit data delete(dux,duy,dvx) finalize
!$acc exit data delete(xpnl,xpne) finalize
!$acc exit data delete(dpnl,dpne) finalize

!$acc exit data copyout(z,xp,dp) finalize

! Export field to a file @ fields/w.* /Xp.* /Dp.*
call writef(z,ifr)

write(nome,"('./fields/Dp.',i3.3)")ifr
open(unit=70,file=nome,form='unformatted')
write(70)dp
close(70)

write(nome,"('./fields/Xp.',i3.3)")ifr
open(unit=70,file=nome,form='unformatted')
write(70)xp
close(70)

! ========================================================

err_inv = err_inv + cufftDestroy(plan_inv)
err_dir = err_dir + cufftDestroy(plan_dir)
istat = cudaStreamDestroy(stream1)
istat = cudaStreamDestroy(stream2)

! Calculate the "mean" time per iteration
! call system_clock(stopTime,trate)
! s=(dble(stopTime-startTime)/dble(dble(npas)*real(trate,8)))
! write(6,*),"Mean time per iteration      = ",s
! call flush(6)

call elapsedtime_s(et_s)
s=dlog(et_s)-dlog(dble(npas))
write(6,*),"Mean time per iteration      = ",dexp(s)
call flush(6)

end