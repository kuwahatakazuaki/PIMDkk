subroutine Init_Velocity
  use Parameters
  use utility, only: gasdev
  implicit none
  real(8) :: vsigma
  !real(8) :: gasdev, vsigma
  integer :: Imode, Iatom
!
!     /*  vsigma: standard deviation of Maxwell distribution  */
!
  do Imode = 1, Nbead
    do Iatom = 1, Natom
      vsigma = dsqrt(1.d0/beta/fictmass(Iatom,Imode))
      vur(1,Iatom,Imode) = vsigma*gasdev()
      vur(2,Iatom,Imode) = vsigma*gasdev()
      vur(3,Iatom,Imode) = vsigma*gasdev()
    enddo
  enddo

  call Remove_TnR_All
  return
end subroutine Init_Velocity
