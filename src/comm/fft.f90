! Direct nd inverse FFT routine
! The direct has the normalization 1/(NX*NY)

subroutine fft_dir(z,plan)
  use cufft
  use openacc
  use paran
  implicit none

  real(8), dimension(NXP2,NY) :: z
  integer :: plan
  integer :: j,ierr
  
  !$acc host_data use_device(z)
  ierr = ierr + cufftExecD2Z(plan,z,z)
  !$acc end host_data
  
  ! In this particular case, single loop 
  ! improves performance because of memory organization
  !$acc parallel loop present(z)
  do j=1,NY
  z(:,j) = z(:,j)/NY/NX
  enddo
  !$acc end parallel loop

return
end

! ========================================================

subroutine fft_inv(z, plan)
  use cufft
  use openacc
  use paran
  implicit none
  
  real(8), dimension(NXP2,NY) :: z
  integer :: plan
  integer :: ierr

  !$acc host_data use_device(z)
  ierr = ierr + cufftExecZ2D(plan,z,z)
  !$acc end host_data

return
end