subroutine inieuler(z,ifr)
  use paran
  use plan
  use commsim
  use commforc
  use commphy
  use commk
  implicit none

  real(8), intent(out   ), dimension(2,NX2P1,NY) :: z
  integer :: ifr,status
  
  open(unit=9,file='./params.dat')
  read(9,*)xlx,xly
  read(9,*)nu,alpnu,mu,alpmu,alpv
  read(9,*)dt,rko
  read(9,*)npas,iout,imix
  read(9,*)rtrunc,truncate
  read(9,*)famp,k1f,k2f
  close(9)
  
  split=0
  
  if (alpnu < 1.0005) then 
  if (alpnu > 0.9995) then 
  	if (alpmu < 0.0005) then 
  	if (alpmu >-0.0005) then 
  		split=1
  	end if
  	end if
  end if
  end if
  
  if (rko==4) then 
  	split=0
	end if
	
  write(6,*)' 2D alpha-turbulence simulation'
  write(6,*)' Resolution       = ',NX,NY
  write(6,*)' Viscosity        = ',real(nu)
  write(6,*)' Viscosity order  = ',real(alpnu)
  write(6,*)' Friction         = ',real(mu)
  write(6,*)' Friction order   = ',real(alpmu)
  write(6,*)' Alpha            = ',real(alpv)
  write(6,*)' '
  
  if (truncate.eq.1) then
   write(6,*)' Truncation at   = ',real(rtrunc)
  else
   write(6,*)' No truncation'
  end if
  
  if (split.eq.1) then
   write(6,*)' Splitted linear integrator (only for Ekman-NSE)'
  else
   write(6,*)' Unsplitted linear integrator'
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
  
  ! Check previous simulations
  if (ifr.gt.0) then
    call readf(z,ifr)
  else
  	! Otherwise, start from zero
    z(:,:,:) = 0.d0
  endif

return
end