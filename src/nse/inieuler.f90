! Code initialization. Print parameter, call initialization
! of wave vector, call initialization of forcing, 
! load / initialize starting field
! Also define log-spaced point to sample structure functions


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
  integer :: iold,is,i,j
  real(8) :: nmax
  
  open(unit=9,file='./params.dat')
  read(9,*)xlx,xly
  read(9,*)nu,alpnu,mu,alpmu,alpv
  read(9,*)dt,rko,prec
  read(9,*)npas,iout,imix
  read(9,*)istru,nskip
  read(9,*)rtrunc,truncate
  read(9,*)famp,k1f,k2f
  close(9)
  	
  write(6,*)' alpha-Turbulence Simulation'
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
  write(6,*)' Out strucs / skip   = ',istru,nskip
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
  
  ! Initialize forcing, see iniforcing.f90
  call iniforcing
  
  ! Check previous simulations
  if (ifr.gt.0) then
    call readf(z,ifr)
  else
  	! Otherwise, start from zero
    z(:,:,:) = 0.d0
  endif
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define log-binned scales !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

iold = 0
is   = 1
j    = 0
nmax = 4 * log(dble(NX / 2)) / log(2.d0)

10 continue
  i = 2**(0.25d0 * dble(is))
  if (i .eq. iold) then
    is = is + 1
    goto 10
  end if
  j = j + 1
  iss(j) = i
  iold = i
  is = is + 1
  if (is .gt. nmax) goto 20
  goto 10
20 nscale = j

return
end