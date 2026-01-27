! Code initialization. Print parameter, call initialization
! of wave vector, call initialization of forcing, 
! load / initialize starting field

subroutine inieuler(z,w,ifr)
  use paran
  use plan
  use commsim
  use commforc
  use commphy
  use commk
  use pass
  implicit none

  real(8), intent(out   ), dimension(2,NX2P1,NY) :: z,w
  integer :: ifr,status
  
  open(unit=9,file='./params.dat')
  read(9,*)xlx,xly
  read(9,*)nu,alpnu,mu,alpmu,alpv
  read(9,*)nup,alpnup,mup,alpmup
  read(9,*)dt,rko,prec
  read(9,*)npas,iout,imix
  read(9,*)rtrunc,truncate
  read(9,*)famp,k1f,k2f
  close(9)
	
  write(6,*)' alphaTurb + Passive scalar'
  write(6,*)' Resolution        = ',NX,NY
  write(6,*)' Viscosities       = ',real(nu),real(nup)
  write(6,*)' Viscosities order = ',real(alpnu),real(alpnup)
  write(6,*)' Frictions         = ',real(mu),real(mup)
  write(6,*)' Frictions order   = ',real(alpmu),real(alpmup)
  write(6,*)' Alpha             = ',real(alpv)
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
  write(6,*)' Saving Precision    = ',prec
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
    call readp(w,ifr)
  else
  	! Otherwise, start from zero
    z(:,:,:) = 0.d0
    w(:,:,:) = 0.d0
  endif

return
end