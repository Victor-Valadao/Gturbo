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