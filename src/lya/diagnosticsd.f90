! Diagnostics for the difference field ~w-z includes:
! spectrum            -> './files/spectradiff.*' (22)
! global dissipations -> './globaldiff.*'        (16)

! ================================================================

subroutine spectrumdiff(z,w,ifr)
use paran
use commk
use commphy
use cudafor
implicit none

complex(8), intent(in), dimension(NX2P1,NY) :: z,w
real(8), dimension(NBIN) :: sp,sz,wn
real(8) k2,ka,k,dk,fac,scra,rho
integer :: i,j,ib,ifr
character(len=64) :: nome

dk=dkx
rho=2.d0/3.14159265358979323846d0
sp=0.d0
sz=0.d0
wn=0.d0     ! measures the number of wavevectors in the shell |k|

!$acc enter data copyin(sp,sz,wn)
!$acc parallel loop collapse(2) present(z,fkxs,fkys,sp,sz,wn)
do 10 j=1,NY
do 10 i=1,NX2P1
  if ((i.eq.1).and.(j.eq.1)) goto 10
  fac=1.d0
  if ((i.eq.1).or.(i.eq.NX2P1)) fac=0.5d0
  k2=fkxs(i) + fkys(j)
  k = sqrt(k2)
  ka=k2**(alpv/2.d0)
  ib = int(k/dk)+1
  scra=0.5d0*cdabs(z(i,j)-w(i,j))**2.0
  if (ib.lt.NBIN) then
 	!$acc atomic
    sp(ib)=sp(ib)+fac*scra/ka
    !$acc atomic
    sz(ib)=sz(ib)+fac*scra
	!$acc atomic
    wn(ib)=wn(ib)+fac*rho/dble(dk*(ib**2-(ib-1)**2))
  end if
10 continue
!$acc end parallel

!$acc kernels present(sp,sz,wn)
sp=sp/dk
sz=sz/dk
wn=wn/dk
wn(1)=1.d0
wn(NBIN)=1.d0
!$acc end kernels

!$acc update host(sp,sz,wn)
do ib=1,NBIN
 write(22,99)real(dk*(ib-1)),real(sp(ib)/wn(ib)),real(sz(ib)/wn(ib))
 99 format(3g)
enddo
 
!$acc exit data delete(sp,sz,wn) finalize

return
end

! ================================================================

subroutine globdiff(z,w,ifr)
use paran
use commk
use commphy
use commforc
implicit none

complex(8), intent(in), dimension(NX2P1,NY) :: z,w
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

!$acc enter data copyin(ee,zz,dise_nu,dise_mu,disz_nu,disz_mu)
!$acc parallel loop collapse(2) present(z,fkxs,fkys,ee,zz,dise_nu,dise_mu,disz_nu,disz_mu)
do 10 j=1,NY
do 10 i=1,NX2P1
 fac=1.d0
 if ((i.eq.1).or.(i.eq.NX2P1)) fac=0.5d0
 if ((i.eq.1).and.(j.eq.1)) goto 10
 k2 = fkxs(i) + fkys(j)
 ka=k2**(alpv/2.d0)
 scra=0.5d0*cdabs(z(i,j)-w(i,j))**2
 !$acc atomic
 ee=ee+fac*scra/ka
 !$acc atomic
 zz=zz+fac*scra
 !$acc atomic
 disz_nu=disz_nu+2.d0*fac*nu*(k2**alpnu)*scra
 !$acc atomic
 disz_mu=disz_mu+2.d0*fac*mu*(k2**alpmu)*scra
 !$acc atomic
 dise_nu=dise_nu+2.d0*fac*nu*(k2**(alpnu-alpv/2.d0))*scra
 !$acc atomic
 dise_mu=dise_mu+2.d0*fac*mu*(k2**(alpmu-alpv/2.d0))*scra
10 continue
!$acc end parallel

!$acc update host(ee,zz,dise_nu,dise_mu,disz_nu,disz_mu)
write(16,99)real(ifr*dt),real(zz),real(ee),real(disz_nu),real(disz_mu),real(inpz),real(dise_nu),real(dise_mu),real(inpe)
call flush(16)
99 format(9g)

!$acc exit data delete(ee,zz,dise_nu,dise_mu,disz_nu,disz_mu) finalize

return
end