subroutine iniforcing
use paran
use commk
use commforc
use commphy
implicit none

real(8) :: k2,k2n
integer :: i,j

nf=0
inpe=0.d0
inpz=0.d0

i=1
do j=2,NY2P1
 k2 = fkxs(i) + fkys(j)
 k2n=fkxs(i)/fkxsmax+fkys(j)/fkysmax
 if (k2n.le.rtrunc) then
 if ((k2.le.k2f**2).and.(k2.gt.k1f**2)) then
  nf=nf+1
  if (nf.gt.NFORC) stop 'NFORC too small !'
  ikxf(nf)=i
  ikyf(nf)=j
  ff(1,nf)=1.d0
  ff(2,nf)=1.d0
  inpz=inpz+famp**2
  inpe=inpe+famp**2/k2**(alpv/2.d0)
 end if
 end if
end do

do i=2,NX2
do j=1,NY
 k2 = fkxs(i) + fkys(j)
 k2n=fkxs(i)/fkxsmax+fkys(j)/fkysmax
 if (k2n.le.rtrunc) then
 if ((k2.le.k2f**2).and.(k2.gt.k1f**2)) then
  nf=nf+1
  if (nf.gt.NFORC) stop 'NFORC too small !'
  ikxf(nf)=i
  ikyf(nf)=j
  ff(1,nf)=1.d0
  ff(2,nf)=1.d0
  inpz=inpz+famp**2
  inpe=inpe+famp**2/k2**(alpv/2.d0)
 end if
 end if
end do
end do

i=NX2P1
do j=1,NY2P1
 k2 = fkxs(i) + fkys(j)
 k2n=fkxs(i)/fkxsmax+fkys(j)/fkysmax
 if (k2n.le.rtrunc) then
 if ((k2.le.k2f**2).and.(k2.gt.k1f**2)) then
  nf=nf+1
  if (nf.gt.NFORC) stop 'NFORC too small !'
  ikxf(nf)=i
  ikyf(nf)=j
  ff(1,nf)=1.d0
  ff(2,nf)=1.d0
  inpz=inpz+famp**2
  inpe=inpe+famp**2/k2**(alpv/2.d0)
 end if
 end if
end do

write(6,*)' # of forced wavenumbers = ',nf
write(6,*)' Mean forcing amplitude  = ',real(famp)
write(6,*)' Forcing range           : ',real(k1f),' < k =<',real(k2f)
write(6,*)' Mean energy input       = ',real(inpe)
write(6,*)' Mean enstrophy input    = ',real(inpz)
write(6,*)' '

call flush(6)

return
end