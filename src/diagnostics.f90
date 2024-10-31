!!!!!!!!!!!!!!!!!!!!!
! Calculate spectra !
!!!!!!!!!!!!!!!!!!!!!

subroutine spectrum(z,ifr)
use paran
use commk
use commphy
use cudafor
implicit none

complex(8), intent(in), dimension(NX2P1,NY) :: z
real(8), dimension(NBIN) :: sp,sz
real*8 k2,ka,k,dk,fac
integer :: i,j,ib,ifr
character*64 :: nome
integer(kind=cuda_count_kind) :: fre,tota

dk=dkx
sp=0.d0
sz=0.d0

!$acc enter data copyin(sp,sz)

!$acc parallel loop collapse(2) present(z,fkxs,fkys,sp,sz)
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
  end if
10 continue
!$acc end parallel

!$acc kernels present(sp,sz)
sp=sp/dk
sz=sz/dk
!$acc end kernels

!$acc update host(sp,sz)
do ib=1,NBIN
 write(20,99)real(dk*(ib-1)),real(sp(ib)),real(sz(ib))
 99 format(3g)
end do

!$acc exit data delete(sp,sz) finalize

return
end

!!!!!!!!!!!!!!!!!!!!!
! Calculate fluxes  !
!!!!!!!!!!!!!!!!!!!!!

subroutine fluxes(z,znl,ifr)
use paran
use commk
use commphy
implicit none

real(8), dimension(2,NX2P1,NY) :: z,znl

real(8), dimension(NBIN) :: fle,flz,flem,flzm
real*8 k2,k,ka,dk,fe,fz,fac
integer i,j,ib,jb,ifr
character*40 nome

dk=dkx

flem=0.d0
flzm=0.d0
fle=0.d0
flz=0.d0

!$acc enter data copyin(flem,flzm,fle,flz)

!$acc parallel loop collapse(2) present(z,znl,fkxs,fkys,flem,flzm)
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
    end if
10 continue
!$acc end parallel

!$acc parallel loop present(fle,flz,flem,flzm)
do i=1,NBIN   ! Flux across the mode k towards large scales
	flz(i)=sum(flzm(i:NBIN))
	fle(i)=sum(flem(i:NBIN))
end do
!$acc end parallel

!$acc update host(flz,fle)
do ib=1,NBIN
 write(30,99)real(dk*(ib)),real(fle(ib)),real(flz(ib))
 99 format(3g)
end do

!$acc exit data delete(flem,flzm,fle,flz) finalize
return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate dissipative terms !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine glob(z,ifr)
use paran
use commk
use commphy
use commforc
implicit none

complex(8), intent(in), dimension(NX2P1,NY) :: z
real(8) :: k2,ka,fac,ee,zz,scra
real(8) :: dise_nu,dise_mu,disz_nu,disz_mu
integer :: i,j,ifr
character*40 :: nome

ee=0.0d0
zz=0.0d0
dise_nu=0.d0
dise_mu=0.d0
disz_nu=0.d0
disz_mu=0.d0

!$acc enter data copyin(ee,zz,dise_nu,dise_mu,disz_nu,disz_mu)

!$acc parallel loop collapse(2) present(z,fkxs,fkys,ee,zz,dise_nu,dise_mu,disz_nu,disz_mu)
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
!$acc end parallel

!$acc update host(ee,zz,dise_nu,dise_mu,disz_nu,disz_mu)
write(18,99)real(ifr*dt),real(zz),real(ee),real(disz_nu),real(disz_mu),real(inpz),real(dise_nu),real(dise_mu),real(inpe)
call flush(18)

99 format(9g)

! write(6,*)"t  = ",real(ifr*dt),", Z  = ",real(zz),", E  = ",real(ee)
write(6,*)"Nt = ",ifr,", Zb = ",real((disz_nu+disz_mu)/inpz),", Eb = ",real((dise_nu+dise_mu)/inpe)
! write(6,*)" "

call flush(6)

!$acc exit data delete(ee,zz,dise_nu,dise_mu,disz_nu,disz_mu) finalize
return
end