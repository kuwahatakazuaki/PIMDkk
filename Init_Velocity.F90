subroutine Init_Velocity
  use Parameters
  use utility, only: fill_gaussian_random
  implicit none
  real(8), allocatable :: sigma(:,:), sigma3(:,:,:)
  !real(8), allocatable :: z(:,:,:)
!
!     /*  sigma: standard deviation of Maxwell distribution  */
!
  allocate(sigma(Natom,Nbead), sigma3(Ndim,Natom,Nbead))
  !allocate(z(Ndim,Natom,Nbead))

  sigma(:,:) = dsqrt(1.d0/beta/fictmass(:,:))
  sigma3(:,:,:) = spread(sigma, dim=1, ncopies=Ndim)
  call fill_gaussian_random(vur)
  vur(:,:,:) = sigma3(:,:,:) * vur(:,:,:)

  call Remove_TnR_All
  deallocate(sigma, sigma3)
  return
end subroutine Init_Velocity
