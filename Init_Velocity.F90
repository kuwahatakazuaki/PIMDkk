subroutine Init_Velocity
  use Parameters
  implicit none
  Double Precision :: gasd, vsigma
  integer :: imode, iatom
!
!     /*  vsigma: standard deviation of Maxwell distribution  */
!
  do imode = 1, nbead
    do iatom = 1, natom
      vsigma = dsqrt(1.d0/beta/fictmass(iatom,imode))
      call gasdev(gasd); vur(1,iatom,imode) = vsigma*gasd
      call gasdev(gasd); vur(2,iatom,imode) = vsigma*gasd
      call gasdev(gasd); vur(3,iatom,imode) = vsigma*gasd
    enddo
  enddo

  call Remove_TnR_All
  return
end subroutine Init_Velocity
