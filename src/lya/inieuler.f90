! Code initialization. Print parameter, call initialization
! of wave vector, call initialization of forcing, 
! load / initialize starting field

subroutine inieuler(z,w,ifr,seed)
  use paran
  use plan
  use commsim
  use commforc
  use commphy
  use commk
  use lyap
  implicit none

  real(8), intent(out   ), dimension(2,NX2P1,NY) :: z,w
  integer :: ifr,status,seed
  
  open(unit=9,file='./params.dat')
  read(9,*)xlx,xly
  read(9,*)nu,alpnu,mu,alpmu,alpv
  read(9,*)dt,rko,prec
  read(9,*)npas,iout,imix,nlyap
  read(9,*)rtrunc,truncate
  read(9,*)famp,k1f,k2f,delta0
  close(9)
	
  write(6,*)' alphaTurb + Eulerian FTLE'
  write(6,*)' Resolution      = ',NX,NY
  write(6,*)' Viscosity       = ',real(nu)
  write(6,*)' Viscosity order = ',real(alpnu)
  write(6,*)' Friction        = ',real(mu)
  write(6,*)' Friction order  = ',real(alpmu)
  write(6,*)' Alpha           = ',real(alpv)
  write(6,*)' '
  
  if (truncate.eq.1) then
   write(6,*)' Truncation at   = ',real(rtrunc)
  else
   write(6,*)' No truncation'
  end if
  
  if (rko.eq.2) then
   write(6,*)' Runge-Kutta of 2nd order'
  else if (rko.eq.4) then
   write(6,*)' Runge-Kutta of 4th order'
   else
   write(6,*)' Runge-Kutta order invalid, falling back to 2nd'
   rko=2
  end if
  write(6,*)' '
  write(6,*)' Time step           = ',real(dt)
  write(6,*)' Starting Frame      = ',ifr
  write(6,*)' Total step number   = ',npas
  write(6,*)' Outputs Frames/Diag = ',iout,imix
  write(6,*)' ------------------------------------------------'
  write(6,*)' '
  
  call flush(6)
	
	! Calculate reusable floats
	sdt=sqrt(dt)
    dt2=dt/2.d0
	dt3=dt/3.d0
	dt6=dt/6.d0
  
  ! Initialize wavenumber, see defk.f90
  call defk
  
  ! Initialize forcing, see iniforcing1.f90
  call iniforcing
  
  if (ifr.eq.0) then ! Generate random noise
!     call readf(z,1)
	z=0.d0
    call add_noise(w,seed)
    w = z + w
    call writep(w,0)
    call writef(z,0)
!     ifr=ifr+1
    
  else   ! Otherwise load previous simulations
    call readf(z,ifr)
    call readp(w,ifr)
  endif

return
end