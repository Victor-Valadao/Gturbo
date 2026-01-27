! Diagnostics for the secondly integrated field, includes:
! spectrum            -> './files/spectrap.*' (21)
! fluxes              -> './files/fluxesp.*'  (31)
! global dissipations -> './globalp.*'        (19)

! ================================================================

subroutine spectrump(z,ifr)
use paran
use commk
use commphy
use cudafor
implicit none

complex(8), intent(in), dimension(NX2P1,NY) :: z
real(8), dimension(NBIN) :: sp,sz
real(8) k2,ka,k,dk,fac
integer :: i,j,ib,ifr
character(len=64) :: nome
integer(kind=cuda_count_kind) :: fre,tota

dk=dkx
sp=0.d0
sz=0.d0

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
    sp(ib)=sp(ib)+fac*k2*abs(z(i,j))**2.0  ! gradient statistics
    sz(ib)=sz(ib)+fac*abs(z(i,j))**2.0
  end if
10 continue

sp=sp/dk
sz=sz/dk

do ib=1,NBIN
 write(21,99)real(dk*(ib-1)),real(sp(ib)),real(sz(ib))
 99 format(3g)
end do

return
end

! ================================================================

subroutine fluxesp(z,znl,ifr)
use paran
use commk
use commphy
implicit none

real(8), dimension(2,NX2P1,NY) :: z,znl

real(8), dimension(NBIN) :: fle,flz,flem,flzm
real(8) k2,k,ka,dk,fe,fz,fac
integer i,j,ib,jb,ifr
character(len=64) :: nome

dk=dkx
fle=0.d0
flz=0.d0
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
    fz=z(1,i,j)*znl(1,i,j)+z(2,i,j)*znl(2,i,j)  ! -T_k / 2
	ka=k2**(alpv/2.d0)
    fe=fz/ka
    ib = int(k/dk)+1
    if (ib.lt.NBIN) then  ! Flux contribution for each mode
    	flzm(ib)=flzm(ib)+fac*fz
    	flem(ib)=flem(ib)+fac*fe
    end if
10 continue

do i=1,NBIN   ! Flux across the mode k towards large scales
	flz(i)=sum(flzm(i:NBIN))
	fle(i)=sum(flem(i:NBIN))
end do

do ib=1,NBIN
 write(31,99)real(dk*(ib)),real(fle(ib)),real(flz(ib))
 99 format(3g)
end do

return
end

! ================================================================

subroutine globp(z,ifr)
use paran
use commk
use commphy
use commforc
implicit none

complex(8), intent(in), dimension(NX2P1,NY) :: z
real(8) :: k2,ka,fac,ee,zz,scra
real(8) :: dise_nu,dise_mu,disz_nu,disz_mu
integer :: i,j,ifr
character(len=64) :: nome

ee=0.0d0
zz=0.0d0

dise_nu=0.d0
dise_mu=0.d0
disz_nu=0.d0
disz_mu=0.d0

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
 disz_mu=disz_mu+2.d0*fac*mu*(k2**alpmu)*scra
 dise_nu=dise_nu+2.d0*fac*nu*(k2**(alpnu-alpv/2.d0))*scra
 dise_mu=dise_mu+2.d0*fac*mu*(k2**(alpmu-alpv/2.d0))*scra
10 continue

write(19,99)real(ifr*dt),real(zz),real(ee),real(disz_nu),real(disz_mu),real(inpz),real(dise_nu),real(dise_mu),real(inpe)
call flush(19)

99 format(9g)

! write(6,*)"t  = ",real(ifr*dt),", Z  = ",real(zz),", E  = ",real(ee)
! write(6,*)"Passive - Nt = ",ifr,", Zb = ",real((disz_nu+disz_mu)/inpz)
! write(6,*)" "

call flush(6)
return
end