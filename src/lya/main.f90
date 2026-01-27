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
use lyap
use time_keeper
implicit none

real(8), dimension(NS) :: z,w,znl,zne,u,v
integer :: ipas,ifr,jfr,istat,seed,seedp,ierr,err_dir,err_inv,i,j
integer :: startTime, stopTime, trate
integer(kind=cuda_count_kind) :: fre,tota
real :: lambda
real(8) :: s,t0,et_s,DeD
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

seedp=seed

! Initializes wave number, forcing, initial field z (real space)
call inieuler(z, w, ifr, seedp)

seedp=seed

! Opening files to save global quantities and useful information
write(nome,"('./files/global.',i3.3)")ifr
Open(unit=18,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/enegspectra.',i3.3)")ifr
Open(unit=14,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/enegdiff.',i3.3)")ifr
Open(unit=15,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/globaldiff.',i3.3)")ifr
Open(unit=16,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/resctime.',i3.3)")ifr
Open(unit=17,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/spectra.',i3.3)")ifr
Open(unit=20,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/fluxes.',i3.3)")ifr
Open(unit=30,file=nome,form='formatted',status='unknown',access='append')
! Opening files to save global quantities and useful information (passive)
write(nome,"('./files/globalp.',i3.3)")ifr
Open(unit=19,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/spectrap.',i3.3)")ifr
Open(unit=21,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/spectradiff.',i3.3)")ifr
Open(unit=22,file=nome,form='formatted',status='unknown',access='append')

write(nome,"('./files/fluxesp.',i3.3)")ifr
Open(unit=31,file=nome,form='formatted',status='unknown',access='append')

! Shows the amount of memory consumed in the GPU
call show_mem(fre,tota,1)

! ========================================================

! Copying data to the GPU data

!$acc enter data copyin(emk)                  !<-- calc in defk called in inieuler
!$acc enter data copyin(ikxf,ikyf,ff)         !<-- calc in iniforcing called in inieuler
!$acc enter data copyin(fkx,fky,fkxs,fkys,jc) !<-- calc in defk called in inieuler
!$acc enter data copyin(dt,dt2,dt3,dt6,sdt)   !<-- calc in inieuler/loaded
!$acc enter data copyin(seed,seedp,eed,delta0)!<-- read in main used in forcing/rann
!$acc enter data copyin(znl,zne,u,v)          !<-- allocate memory for znl, u, v

!$acc enter data copyin(z,w)                  !<-- pass the initial field to the GPU

! Shows the amount of memory consumed in the GPU
call show_mem(fre,tota,2)

if (ifr.eq.0) then ! Generate random noise
	call fft_dir(z,plan_dir)
	call perturb(z,w,eta*xlx/NX,beta,seed)
	if (truncate.eq.1) then
	 call trunc(w)
	end if
	call fft_inv(w,plan_inv)
	call fft_inv(z,plan_inv)
	!$acc update host(w,z)
	call writep(w,0)
end if

! DirFFT the initial field + dealiasing
call fft_dir(z,plan_dir)
call fft_dir(w,plan_dir)
if (truncate.eq.1) then
 call trunc(z)
 call trunc(w)
end if

DeD=0.d0
! $acc update host(z,w)
call enediff(z,w,DeD,0)
call enespec(z,w,0)
! call rescale(z,w,eed)
call spectrumdiff(z,w,0)
! call globdiff(z,w,0)

jfr=1
! =============================================================
!  Time evolution <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
do ipas=1,npas
  t=t+dt
  
  if (rko.eq.4) then
  	call step4(z,znl,zne,u,v)
  else
  	call step(z,znl,u,v)
  end if
  
  if (rko.eq.4) then
  	call step4(w,znl,zne,u,v)
  else
  	call step(w,znl,u,v)
  end if
  
  if (mod(ipas,imix).eq.0) then  ! Write diagnostics each imix
    !$acc kernels present(z,w,znl,zne) ! Copy z to znl
	znl = z
	zne = w
	!$acc end kernels
	call nlt(znl ,u,v)            ! nlt to compute the flux at t correctly
	call nlt(zne ,u,v)
		
  	! $acc update host(z,znl,w,zne)
  	call spectrum(z,ipas)        ! Write Energy and Enstrophy spectra
    call fluxes(z,znl,ipas)      ! Write Energy and Enstrophy fluxes
    call glob(z,ipas)            ! Write Dissipative terms
	call spectrump(w,ipas)       ! Write Energy and Enstrophy spectra
    call fluxesp(w,zne,ipas)     ! Write Energy and Enstrophy fluxes
    call globp(w,ipas)
	call spectrumdiff(z,w,ipas)
    call globdiff(z,w,ipas)
  endif
  
  if ((mod(ipas,iout).eq.0).and.(ipas.ne.npas)) then  ! Save vorticity each iout
    !$acc kernels present(z,w,znl,zne) ! Copy z to znl
	znl = z
	zne = w
	!$acc end kernels
	
	call fft_inv(znl , plan_inv)  ! FFt to save in the physical space
	call fft_inv(zne , plan_inv)
	!$acc update host(znl,zne) 
  	
  	call writef2(znl ,ifr,jfr)        ! Save field
  	call writep2(zne ,ifr,jfr)
  	jfr=jfr+1
  	
  	write(6,*)"Saved w field, seed and curframe at timestep = ",ipas
  	call flush(6)
  	
  endif
  
	! ----------------------------------------------------------
	! Lyapunov
	if (mod(ipas,nlyap).eq.0) then 
	 Ded=0.d0
	 !$acc kernels present(z,w,znl,zne) ! Copy z to znl
		znl = z
		zne = w
	 !$acc end kernels
	 
	 call enediff(znl,zne,DeD,ipas)
!      call enespec(z,w,0)
	 
	 if (DeD.gt.4) then
! 	 	 write(6,*)"Vrau cond on em neguin"
! 		 call flush(6)
		 		 
		 call enespec(znl,zne,ipas)
		 write(17,*)ipas
		 call perturb(znl,zne,eta*xlx/NX,beta,seed)
		write(6,*)"RESCALED EM KRL PORRA EM ",ipas
		call flush(6)
		 
		 if (truncate.eq.1) then
		 call trunc(znl)
		 call trunc(zne)
		 
		!$acc kernels present(zne,w)
		w = zne
		!$acc end kernels
! 	
! 		call fft_inv(znl , plan_inv)  ! FFt to save in the physical space
! 		call fft_inv(zne , plan_inv)
! 		!$acc update host(znl,zne) 
! 	
! 		call writef2(znl ,ifr,jfr)        ! Save field
! 		call writep2(zne ,ifr,jfr)
! 		jfr=jfr+1
! 	
! 		write(6,*)"Saved w field, seed and curframe at timestep = ",ipas
! 		call flush(6)
		 
		end if
		 
	  end if
	end if
	! ----------------------------------------------------------
  
    ! ! Add forcing contribution
    !$acc kernels present(seed,seedp)
	seedp=seed
	!$acc end kernels  
  call forcing(z,seed)
  call forcing(w,seedp)
   
  
end do
! --------------------------------------

ifr=ifr+1
! Final write
!$acc update host(z,w)
call fft_inv(z, plan_inv)
call fft_inv(w, plan_inv)

!$acc update host(seed,seedp)
write(nome,"('./files/seed.',i3.3)")ifr
open(unit=2,file=nome)
write(2,*)seed
close(2)

! Export data to curframe.dat to continue, if you want
open(unit=1,file='./curframe.dat')
 write(1,*)ifr,t
close(1)

close(14)
close(15)
close(16)
close(17)
close(18)
close(20)
close(30)
close(19)
close(21)
close(22)
close(31)

! ========================================================

! Exit GPU data Freeing the memory

!$acc exit data delete(emk) finalize
!$acc exit data delete(ikxf,ikyf,ff) finalize
!$acc exit data delete(fkx,fky,fkxs,fkys,jc) finalize
!$acc exit data delete(dt,dt2,dt3,dt6,sdt) finalize
!$acc exit data delete(seed,seedp,eed,delta0) finalize
!$acc exit data delete(u,v,znl,zne) finalize

!$acc exit data copyout(z,w) finalize

! Export field to a file @ Frames/w.*
call writef(z,ifr)
call writep(w,ifr)
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