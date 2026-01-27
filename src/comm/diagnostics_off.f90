! Diagnostics for the main integrated field includes:
! spectrum            -> './files/spectra.*' (20)
! fluxes              -> './files/fluxes.*'  (30)
! global dissipations -> './global.*'        (18)

! ================================================================

subroutine spectrum(z,ifr)
use paran
use commk
use commphy
use cudafor
implicit none

complex(8), intent(in), dimension(NX2P1,NY) :: z
real(8), dimension(NBIN) :: sp,sz,wn
real(8) k2,ka,k,dk,fac,rho
integer :: i,j,ib,ifr
character(len=64) :: nome
integer(kind=cuda_count_kind) :: fre,tota

dk=dkx
rho=2.d0/3.14159265358979323846d0
sp=0.d0
sz=0.d0
wn=0.d0     ! measures the number of wavevectors in the shell |k|

do 10 j=1,NY
do 10 i=1,NX2P1
  if ((i.eq.1).and.(j.eq.1)) goto 10
  fac=1.d0
  if ((i.eq.1).or.(i.eq.NX2P1)) fac=0.5d0
  k2=fkxs(i) + fkys(j)
  k = sqrt(k2)
  ka=k2**(alpv/2.d0)
  ib = int(k/dk)+1
  if (ib.lt.NBIN) then
    sp(ib)=sp(ib)+fac*abs(z(i,j))**2.0/ka
    sz(ib)=sz(ib)+fac*abs(z(i,j))**2.0
    wn(ib)=wn(ib)+fac*rho/dble(dk*(ib**2-(ib-1)**2))
  end if
10 continue

sp=sp/dk
sz=sz/dk
wn=wn/dk
wn(1)=1.d0
wn(NBIN)=1.d0

do ib=1,NBIN
 write(20,99)real(dk*(ib-1)),real(sp(ib)/wn(ib)),real(sz(ib)/wn(ib))
 99 format(3g)
!  prints also wn's
!  write(20,99)real(dk*(ib-1)),real(sp(ib)/wn(ib)),real(sz(ib)/wn(ib)),real(wn(ib))
!  99 format(4g)
end do

return
end

! ================================================================

subroutine fluxes(z,znl,ifr)
use paran
use commk
use commphy
implicit none

real(8), dimension(2,NX2P1,NY) :: z,znl

real(8), dimension(NBIN) :: fle,flz,flem,flzm,wn,wle
real(8) k2,k,ka,dk,fe,fz,fac,rho
integer i,j,ib,jb,ifr
character(len=64) :: nome

dk=dkx
rho=2.0d0/3.14159265358979323846d0
fle=0.d0
flz=0.d0
wle=0.d0

wn=0.d0
flem=0.d0
flzm=0.d0

do 10 j=1,NY
  do 10 i=1,NX2P1
!   if ((mod(j,100).eq.0).and.(mod(i,NX2P1).eq.0)) print*,j
    if ((i.eq.1).and.(j.eq.1)) goto 10
    fac=2.d0
    if ((i.eq.1).or.(i.eq.NX2P1)) fac=1.d0
    k2=fkxs(i)+fkys(j)
    k=sqrt(k2)
    fz=z(1,i,j)*znl(1,i,j)+z(2,i,j)*znl(2,i,j)
    ka=k2**(alpv/2.d0)
    fe=fz/ka
    ib = int(k/dk)+1
    if (ib.lt.NBIN) then  ! Flux contribution for each mode
    	flzm(ib)=flzm(ib)+fac*fz
    	flem(ib)=flem(ib)+fac*fe
    	wn(ib)=wn(ib)+fac/2.d0
    end if
10 continue

do i=1,NBIN   ! Flux across the mode k towards large scales
	flz(i)=sum(flzm(i:NBIN))
	fle(i)=sum(flem(i:NBIN))
	wle(i)=(1+sum(wn(1:i)))/(dk*dble(i-1))**2
end do

wle(1)=1.d0
wle(NBIN)=1.d0

do ib=1,NBIN
 write(30,99)real(dk*(ib-1)),real(fle(ib)/wle(ib)/rho),real(flz(ib)/wle(ib)/rho)
 99 format(3g)
end do

return
end

! ================================================================

subroutine glob(z,ifr)
use paran
use commk
use commphy
use commforc
implicit none

complex(8), intent(in), dimension(NX2P1,NY) :: z
real(8) :: k2,ka,fac,ee,zz,scra,km2
real(8) :: dise_nu,dise_mu,disz_nu,disz_mu
integer :: i,j,ifr
character(len=64) :: nome

ee=0.0d0
zz=0.0d0

dise_nu=0.d0
dise_mu=0.d0
disz_nu=0.d0
disz_mu=0.d0

! km2=(k1f+k2f)/2.d0
km2=k1f

do 10 j=1,NY
do 10 i=1,NX2P1
 fac=1.d0
 if ((i.eq.1).or.(i.eq.NX2P1)) fac=0.5d0
 if ((i.eq.1).and.(j.eq.1)) goto 10
 k2 = fkxs(i) + fkys(j)
 ka=k2**(alpv/2.d0)
 scra=cdabs(z(i,j))**2
 ee=ee+fac*scra/ka
 zz=zz+fac*scra
 disz_nu=disz_nu+2.d0*fac*nu*(k2**alpnu)*scra
 dise_nu=dise_nu+2.d0*fac*nu*(k2**(alpnu-alpv/2.d0))*scra
 disz_mu=disz_mu+2.d0*fac*mu*(k2**alpmu)*scra
 dise_mu=dise_mu+2.d0*fac*mu*(k2**(alpmu-alpv/2.d0))*scra
!  if (sqrt(k2).le.km2) then
!    disz_mu=disz_mu+2.d0*fac*mu*(k2**alpmu)*scra
!    dise_mu=dise_mu+2.d0*fac*mu*(k2**(alpmu-alpv/2.d0))*scra
!  end if
   
10 continue

write(18,99)real(ifr*dt),real(zz),real(ee),real(disz_nu),real(disz_mu),real(inpz),real(dise_nu),real(dise_mu),real(inpe)
call flush(18)

99 format(9g)

! write(6,*)"t  = ",real(ifr*dt),", Z  = ",real(zz),", E  = ",real(ee)
write(6,*)"Nt = ",ifr,", Zb = ",real((disz_nu+disz_mu)/inpz),", Eb = ",real((dise_nu+dise_mu)/inpe)
! write(6,*)" "

call flush(6)
return
end