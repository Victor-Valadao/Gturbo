! Module with all parameters of the simulation

! Sim parameters
module paran
 implicit none
 integer, parameter :: NX=8192, NY=NX
 integer, parameter :: NXP2=NX+2
 integer, parameter :: NX2=NX/2, NX2P1=NX2+1, NX2P2=NX2+2
 integer, parameter :: NY2=NY/2, NY2P1=NY2+1, NY2P2=NY2+2
 integer, parameter :: NS=NXP2*NY
 integer, parameter :: NBIN=NX/2
end module paran

! FFT parameters and plans
module stream
  use cufft
  use cudafor
  implicit none
  
  integer(kind=cuda_stream_kind) :: stream1
  integer(kind=cuda_stream_kind) :: stream2
end module stream

module plan
	implicit none

	integer plan_dir
	integer plan_inv
end module plan

! physical parameters 
module commphy
	use paran
	implicit none
	real(8) :: xlx,xly,dkx,dky,prec
	real(8) :: nu,mu,alpnu,alpmu,alpv
	real(8) :: dt,dt2,dt3,dt6,sdt,t
	real(8) :: inpz,inpe,nscale
	integer, dimension(NX) :: iss
end module commphy

! forcing parameters
module commforc
	implicit none
	integer, parameter :: NFORC = 2048
	real(8), dimension(2,NFORC) :: ff
	real(8) :: famp,k1f,k2f
	integer, dimension(NFORC) :: ikxf, ikyf
	integer :: nf
end module commforc

! diagnostics parameters
module commsim
	implicit none
	integer :: npas,iout,imix,rko,oitrp
	integer :: istru,nskip
end module commsim

! wavenumber parameters
module commk
	use paran
	implicit none
	real(8), dimension(NX2P1,NY) :: emk,epk
	real(8), dimension(NX2P1) :: fkxs, fkx
	real(8), dimension(NY) :: fkys, fky
	integer, dimension(NY) :: jc
	real(8) :: rtrunc,fkxsmax,fkysmax
	integer :: truncate
end module commk

! lyap parameters
module lyap
	implicit none
	real(8) :: delta0,eed
	integer :: nlyap
end module lyap

! passive parameters
module pass
	implicit none
	real(8) :: nup,mup,alpnup,alpmup
end module pass

! Lagragian parameters
module modlag
implicit none
	integer, parameter :: NP=100,NV=2
	real(8) :: scalax,scalay,eed
	real(8) :: lyap(NV),ftle(NP,NV)
	real(8) :: co(2,3,3)
	integer :: nlyal
end module

module time_keeper
implicit none
	integer :: start(8), now(8)
contains
  subroutine startclock( )
    implicit none
    call date_and_time(values=start)
  end subroutine startclock

  subroutine elapsedtime_s( et_s )
    implicit none
    integer  :: diffs(8)=0
    real(8)   , intent(out):: et_s       ! in seconds
    call date_and_time(values=now)
    ! - Find the difference in times
    diffs = now - start
    ! - This works only when the time is measured in a specific month
    if (diffs(3) > 0) then
       diffs(5) = 24*diffs(3) + diffs(5)
    endif
    et_s = diffs(5) * 3600 + diffs(6) * 60 + diffs(7) + 1e-3 * diffs(8)
  end subroutine elapsedtime_s
end module time_keeper
