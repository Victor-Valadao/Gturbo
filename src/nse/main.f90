! Main computation loop

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
use time_keeper
implicit none

real(8), dimension(NS) :: z,znl,znew,u,v
integer :: ipas,ifr,jfr,istat,seed,ierr,err_dir,err_inv,i,j
integer :: startTime, stopTime, trate
integer(kind=cuda_count_kind) :: fre,tota
real(8) :: s,t0,et_s
character(len=40) nome

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

! Initializes wave number, forcing, initial field z (real space)
call inieuler(z, ifr)

! Opening files to save global quantities and useful information
write(nome,"('./files/global.',i3.3)")ifr
Open(unit=18,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/spectra.',i3.3)")ifr
Open(unit=20,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/fluxes.',i3.3)")ifr
Open(unit=30,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/strucs.',i3.3)")ifr
Open(unit=40,file=nome,form='formatted',status='unknown',access='append')

! Shows the amount of memory consumed in the GPU
call show_mem(fre,tota,1)

! ========================================================

! Copying data to the GPU data

!$acc enter data copyin(emk,iss)                  !<-- calc in defk called in inieuler
!$acc enter data copyin(ikxf,ikyf,ff)         !<-- calc in iniforcing called in inieuler
!$acc enter data copyin(fkx,fky,fkxs,fkys,jc) !<-- calc in defk called in inieuler
!$acc enter data copyin(dt,dt2,dt3,dt6,sdt)   !<-- calc in inieuler/loaded
!$acc enter data copyin(seed)                 !<-- read in main used in forcing/rann
!$acc enter data copyin(znl,u,v)              !<-- allocate memory for znl, u, v

if (rko.eq.4) then
	!$acc enter data copyin(znew)
endif

!$acc enter data copyin(z)                    !<-- pass the initial field to the GPU

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
  	call step4(z,znl,znew,u,v)
  else
  	call step(z,znl,u,v)
  end if
  
  if (mod(ipas,imix).eq.0) then  ! Write diagnostics each imix
    !$acc kernels present(z,znl) ! Copy z to znl
	znl=z
	!$acc end kernels
	call nlt(znl,u,v)            ! nlt to compute the flux at t correctly
		
  	call spectrum(z,ipas)        ! Write Energy and Enstrophy spectra
    call fluxes(z,znl,ipas)      ! Write Energy and Enstrophy fluxes
    call glob(z,ipas)            ! Write Dissipative terms
	
  endif	
  if (mod(ipas,istru).eq.0) then
	!$acc kernels present(z,znl) ! Copy z to znl
	znl=z
	!$acc end kernels
	call vel(znl,u,v)
	call fft_inv(znl, plan_inv)
	call fft_inv(u, plan_inv)
	call fft_inv(v, plan_inv)
	
	call strucs(znl,u,v)
	
  endif
  
  if ((mod(ipas,iout).eq.0).and.(ipas.ne.npas)) then  ! Save vorticity each iout
  	!$acc kernels present(z,znl) ! Copy z to znl to save without touching z
	znl=z
	!$acc end kernels  
		
	call fft_inv(znl, plan_inv)  ! FFt to save in the physical space
	!$acc update host(znl,seed) 
  	
  	call writef2(znl,ifr,jfr)        ! Save field
  	jfr=jfr+1
  	
  	write(6,*)"Saved w field, seed and curframe at timestep = ",ipas
  	call flush(6)
  	
  endif
  
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

close(18)
close(20)
close(30)
close(40)

! ========================================================

! Exit GPU data Freeing the memory

!$acc exit data delete(emk,iss) finalize
if (rko.eq.4) then
	!$acc exit data delete(znew) finalize
endif

!$acc exit data delete(ikxf,ikyf,ff) finalize
!$acc exit data delete(fkx,fky,fkxs,fkys,jc) finalize
!$acc exit data delete(dt,dt2,dt3,dt6,sdt) finalize
!$acc exit data delete(seed) finalize
!$acc exit data delete(u,v,znl) finalize

!$acc exit data copyout(z) finalize

! Export field to a file @ Frames/w.*
call writef(z,ifr)
! ========================================================

err_inv = err_inv + cufftDestroy(plan_inv)
err_dir = err_dir + cufftDestroy(plan_dir)
istat = cudaStreamDestroy(stream1)
istat = cudaStreamDestroy(stream2)

call elapsedtime_s(et_s)
s=dlog(et_s)-dlog(dble(npas))
write(6,*),"Mean time per iteration      = ",dexp(s)
call flush(6)

end